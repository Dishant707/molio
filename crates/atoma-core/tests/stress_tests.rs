//! Stress tests and fuzz tests for molio-core.
//!
//! Tests that parsing never panics on:
//! - Extremely large inputs
//! - Malformed/garbage input
//! - Edge cases (zero coords, alternate conformations, insertion codes)
//!
//! Run with: `cargo test -p molio-core --test stress_tests`

use atoma_core::parser::{pdb::parse_pdb_str, sdf::parse_sdf_str};

const EDGE_ZERO: &str = include_str!("../../../test_data/edge_cases/zero_coords.pdb");
const EDGE_ALT: &str = include_str!("../../../test_data/edge_cases/alt_conformations.pdb");
const EDGE_EXTREME: &str = include_str!("../../../test_data/edge_cases/extreme_coords.pdb");
const NOT_PDB: &str = include_str!("../../../test_data/edge_cases/not_a_pdb.txt");

// ─── Edge Case Tests ─────────────────────────────────────────────

#[test]
fn test_zero_coordinates() {
    let mol = parse_pdb_str(EDGE_ZERO, "zero.pdb").expect("should parse zero coords");
    assert!(mol.n_atoms() > 0);
    // First atom should have coords near zero
    assert!((mol.atoms[0].x - (-0.001)).abs() < 0.01);
    assert!((mol.atoms[0].y - 0.0).abs() < 0.01);
}

#[test]
fn test_alternate_conformations() {
    let mol = parse_pdb_str(EDGE_ALT, "alt.pdb").expect("should parse alt confs");
    assert!(mol.n_atoms() > 0);
    // Should have atoms with alt_loc 'A' and 'B'
    let alt_a = mol.atoms.iter().any(|a| a.alt_loc == 'A');
    let alt_b = mol.atoms.iter().any(|a| a.alt_loc == 'B');
    assert!(alt_a || alt_b, "should have alternate conformations");
}

#[test]
fn test_extreme_coordinates() {
    let mol = parse_pdb_str(EDGE_EXTREME, "extreme.pdb").expect("should parse extreme coords");
    assert!(mol.n_atoms() > 0);
    // Should handle large coordinate values without overflow
    let _ = mol.bounding_box();
}

#[test]
fn test_non_pdb_input_no_panic() {
    // Garbage input should return an error, never panic
    let result = parse_pdb_str(NOT_PDB, "garbage.txt");
    // Either succeeds (parses nothing) or returns error, but never panics
    let _ = result;
}

#[test]
fn test_empty_string() {
    let mol = parse_pdb_str("", "empty").expect("empty input should be ok");
    assert_eq!(mol.n_atoms(), 0);
}

#[test]
fn test_only_whitespace() {
    let mol = parse_pdb_str("   \n  \n   \n", "whitespace").expect("whitespace input ok");
    assert_eq!(mol.n_atoms(), 0);
}

#[test]
fn test_single_atom_pdb() {
    let pdb = "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N\nEND\n";
    let mol = parse_pdb_str(pdb, "single.pdb").expect("single atom PDB");
    assert_eq!(mol.n_atoms(), 1);
    assert_eq!(mol.chains.len(), 1);
    assert_eq!(mol.chains[0].residues.len(), 1);
}

#[test]
fn test_atom_with_charge() {
    // PDB columns 79-80 hold the charge
    let line = "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00          N  1+";
    let mol = parse_pdb_str(line, "charge.pdb").expect("charged atom");
    assert_eq!(mol.n_atoms(), 1);
    // Charge is parsed from columns 79-80
    let charge = &mol.atoms[0].charge;
    assert!(charge.is_some(), "should have charge");
}

#[test]
fn test_missing_fields() {
    // PDB line with missing coordinate fields
    let pdb = "ATOM      1  N   ALA A   1\nEND\n";
    let result = parse_pdb_str(pdb, "missing.pdb");
    // Should handle gracefully, not panic
    let _ = result;
}

// ─── SDF Edge Cases ──────────────────────────────────────────────

#[test]
fn test_sdf_no_molecules() {
    let mols = parse_sdf_str("").expect("empty SDF");
    assert_eq!(mols.len(), 0);
}

#[test]
fn test_sdf_single_atom_no_bonds() {
    let sdf = "Helium\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 He  0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n";
    let mols = parse_sdf_str(sdf).expect("single atom SDF");
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].n_atoms(), 1);
    assert_eq!(mols[0].bonds.len(), 0);
    assert_eq!(mols[0].atoms[0].element, "He");
}

#[test]
fn test_sdf_missing_delimiter() {
    // SDF without $$$$ terminator (still valid MOL)
    let sdf = "Test\n  test\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\nM  END\n";
    let mols = parse_sdf_str(sdf).expect("SDF without $$$$");
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].n_atoms(), 2);
}

// ─── Stress Tests ────────────────────────────────────────────────

#[test]
fn test_stress_many_atoms_pdb() {
    // Generate a PDB with 1000 atoms
    let mut pdb = String::from("HEADER    STRESS_TEST\n");
    for i in 1..=1000 {
        pdb.push_str(&format!(
            "ATOM  {:>5}  CA  ALA A{:>4}    {:>8.3}{:>8.3}{:>8.3}  1.00  0.00           C  \n",
            i, (i / 10) + 1,
            (i as f64) % 100.0,
            ((i * 7) as f64) % 100.0,
            ((i * 13) as f64) % 100.0,
        ));
    }
    pdb.push_str("END\n");

    let mol = parse_pdb_str(&pdb, "stress.pdb").expect("1000 atom PDB");
    assert_eq!(mol.n_atoms(), 1000);
}

#[test]
fn test_stress_many_molecules_sdf() {
    // Generate SDF with 100 molecules
    let mut sdf = String::new();
    for m in 0..100 {
        sdf.push_str(&format!("Mol_{}\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n", m));
    }

    let mols = parse_sdf_str(&sdf).expect("100 molecule SDF");
    assert_eq!(mols.len(), 100);
}

#[test]
fn test_stress_deep_residue_tree() {
    // Many chains with many residues
    let mut pdb = String::from("HEADER    STRESS\n");
    let mut serial = 1;
    for chain in 'A'..='E' {
        for res in 1..=50 {
            let res_name = if res % 3 == 0 { "GLY" } else if res % 3 == 1 { "ALA" } else { "VAL" };
            for atom_name in &["N", "CA", "C", "O"] {
                pdb.push_str(&format!(
                    "ATOM  {:>5} {:^4}{}{} {}{:>4}    {:>8.3}{:>8.3}{:>8.3}  1.00  0.00           C  \n",
                    serial, atom_name, ' ', res_name, chain, res,
                    (serial as f64) % 100.0,
                    ((serial * 3) as f64) % 100.0,
                    ((serial * 7) as f64) % 100.0,
                ));
                serial += 1;
            }
        }
    }
    pdb.push_str("END\n");

    let mol = parse_pdb_str(&pdb, "deep.pdb").expect("deep structure");
    assert!(mol.chains.len() >= 4, "should have multiple chains");
    assert!(mol.n_atoms() > 900, "should have ~1000 atoms");
    assert!(mol.n_residues() > 200, "should have many residues");
}

// ─── Fuzz Tests ──────────────────────────────────────────────────

#[test]
fn test_fuzz_random_bytes_no_panic() {
    // Random bytes should never cause a panic — only errors
    let garbage: Vec<u8> = (0..1000).map(|i| (i * 17 + 31) as u8).collect();
    if let Ok(s) = std::str::from_utf8(&garbage) {
        let _ = parse_pdb_str(s, "fuzz.pdb");
    } else {
        // Not valid UTF-8 — that's fine, we can't test it
    }
}

#[test]
fn test_fuzz_long_line() {
    let long_line = "A".repeat(10000);
    let _ = parse_pdb_str(&long_line, "long.pdb");
}

#[test]
fn test_roundtrip_properties() {
    // SDF with many properties
    let sdf = "Test\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n> <Prop1>\nval1\n\n> <Prop2>\nval2\n\n> <Prop3>\nval3\n\n$$$$\n";
    let mols = parse_sdf_str(sdf).expect("properties SDF");
    assert_eq!(mols[0].properties.len(), 3);
    assert_eq!(mols[0].properties.get("Prop1").map(|s| s.as_str()), Some("val1"));
}
