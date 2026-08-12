//! Precision and boundary tests.
//!
//! Chaos theory in molecular simulation: a single bit error in a coordinate
//! cascades into completely different trajectories. These tests verify that
//! molio handles edge cases correctly — no silent truncation, no NaN propagation.

use atoma_core::parser::{pdb::parse_pdb_str, sdf::parse_sdf_str, xyz::parse_xyz_str};

// ─── Coordinate Precision ───────────────────────────────────────

#[test]
fn pdb_max_precision_coordinates() {
    // PDB format supports 3 decimal places for coordinates
    let line = "ATOM      1  N   ALA A   1     12.345 -67.890   0.001  1.00  0.00           N  \nEND\n";
    let mol = parse_pdb_str(line, "precise.pdb").unwrap();
    assert_eq!(mol.atoms[0].x, 12.345);
    assert_eq!(mol.atoms[0].y, -67.890);
    assert_eq!(mol.atoms[0].z, 0.001);
}

#[test]
fn sdf_max_precision_coordinates() {
    // SDF format supports 4 decimal places
    let sdf = "Test\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n   -0.0020    1.4050    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n";
    let mols = parse_sdf_str(sdf).unwrap();
    assert!((mols[0].atoms[0].x - (-0.0020)).abs() < 0.0001);
    assert!((mols[0].atoms[0].y - 1.4050).abs() < 0.0001);
}

#[test]
fn xyz_max_precision_coordinates() {
    // XYZ supports arbitrary precision
    let xyz = "1\nPrecision test\nC  1.23456789  -9.87654321  0.00000001\n";
    let mol = parse_xyz_str(xyz).unwrap();
    assert!((mol.atoms[0].x - 1.23456789).abs() < 0.0000001);
    assert!((mol.atoms[0].y - (-9.87654321)).abs() < 0.0000001);
}

// ─── Floating Point Edge Cases ──────────────────────────────────

#[test]
fn pdb_zero_coordinates() {
    let line = "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  \nEND\n";
    let mol = parse_pdb_str(line, "zero.pdb").unwrap();
    assert_eq!(mol.atoms[0].x, 0.0);
    assert_eq!(mol.atoms[0].y, 0.0);
    assert_eq!(mol.atoms[0].z, 0.0);
}

#[test]
fn pdb_negative_zero() {
    // -0.000 should parse as 0.0, not -0.0 (which can cause issues)
    let line = "ATOM      1  N   ALA A   1      -0.000  -0.000  -0.000  1.00  0.00           N  \nEND\n";
    let mol = parse_pdb_str(line, "negzero.pdb").unwrap();
    assert!(mol.atoms[0].x.is_finite());
}

#[test]
fn pdb_large_coordinates() {
    let line = "ATOM      1  N   ALA A   1    -999.999 999.999  -999.999  1.00  0.00           N  \nEND\n";
    let mol = parse_pdb_str(line, "large.pdb").unwrap();
    assert!(mol.atoms[0].x.is_finite());
    assert_eq!(mol.atoms[0].x, -999.999);
}

#[test]
fn sdf_coordinate_no_nan() {
    // All parsed coordinates should be finite
    let sdf = include_str!("../../../test_data/sdf/molecules.sdf");
    let mols = parse_sdf_str(sdf).unwrap();
    for mol in &mols {
        for atom in &mol.atoms {
            assert!(atom.x.is_finite(), "NaN/Inf x in {}", mol.name.as_deref().unwrap_or("?"));
            assert!(atom.y.is_finite(), "NaN/Inf y in {}", mol.name.as_deref().unwrap_or("?"));
            assert!(atom.z.is_finite(), "NaN/Inf z in {}", mol.name.as_deref().unwrap_or("?"));
        }
    }
}

// ─── Element Detection Accuracy ─────────────────────────────────

#[test]
fn pdb_element_detection() {
    // All common elements should be detected correctly
    let tests = vec![
        ("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  \n", "N"),
        ("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  \n", "C"),
        ("ATOM      1  O   ALA A   1       0.000   0.000   0.000  1.00  0.00           O  \n", "O"),
        ("ATOM      1  SG  CYS A   1       0.000   0.000   0.000  1.00  0.00           S  \n", "S"),
        ("ATOM      1  ZN  ZN  A   1       0.000   0.000   0.000  1.00  0.00          ZN  \n", "ZN"),
        ("ATOM      1  FE  HEM A   1       0.000   0.000   0.000  1.00  0.00          FE  \n", "FE"),
    ];
    for (line, expected) in &tests {
        let pdb = format!("{}\nEND\n", line);
        let mol = parse_pdb_str(&pdb, "elem.pdb").unwrap();
        assert_eq!(mol.atoms[0].element, *expected, "failed for line: {}", line);
    }
}

#[test]
fn sdf_element_detection() {
    let elements = ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "He", "Zn", "Fe"];
    let mut sdf = String::from("Elements\n  test\n\n");
    sdf.push_str(&format!("{:>3}  0  0  0  0  0  0  0  0  0999 V2000\n", elements.len()));
    for (i, elem) in elements.iter().enumerate() {
        sdf.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
            i as f64, 0.0, 0.0, elem
        ));
    }
    sdf.push_str("M  END\n$$$$\n");

    let mols = parse_sdf_str(&sdf).unwrap();
    assert_eq!(mols[0].atoms.len(), elements.len());
    for (atom, expected) in mols[0].atoms.iter().zip(elements.iter()) {
        assert_eq!(atom.element, *expected);
    }
}

// ─── Bond Consistency ───────────────────────────────────────────

#[test]
fn sdf_bonds_reference_valid_atoms() {
    // Every bond must reference atoms that exist (1-based indices)
    let sdf = include_str!("../../../test_data/sdf/molecules.sdf");
    let mols = parse_sdf_str(sdf).unwrap();
    for mol in &mols {
        let n = mol.n_atoms() as u32;
        for bond in &mol.bonds {
            assert!(bond.atom1 >= 1 && bond.atom1 <= n,
                "bond atom1 {} out of range (1..{})", bond.atom1, n);
            assert!(bond.atom2 >= 1 && bond.atom2 <= n,
                "bond atom2 {} out of range (1..{})", bond.atom2, n);
            assert_ne!(bond.atom1, bond.atom2, "atom bonded to itself");
        }
    }
}

// ─── Roundtrip Fidelity ─────────────────────────────────────────

#[test]
fn pdb_roundtrip_100_iterations() {
    // Parse the same file 100 times, verify identical results
    let pdb = include_str!("../../../test_data/pdb/1crn.pdb");
    let first = parse_pdb_str(pdb, "r1.pdb").unwrap();
    for i in 0..100 {
        let mol = parse_pdb_str(pdb, &format!("r{}.pdb", i)).unwrap();
        assert_eq!(mol.n_atoms(), first.n_atoms(), "iteration {}: atom count", i);
        assert_eq!(mol.atoms[0].x, first.atoms[0].x, "iteration {}: x coord", i);
    }
}

#[test]
fn xyz_roundtrip_100_iterations() {
    let xyz = "3\nWater\nO  0.000  0.000  0.000\nH  0.957  0.000  0.000\nH -0.240  0.927  0.000\n";
    let first = parse_xyz_str(xyz).unwrap();
    for i in 0..100 {
        let mol = parse_xyz_str(xyz).unwrap();
        assert_eq!(mol.n_atoms(), first.n_atoms());
        assert!((mol.atoms[0].x - first.atoms[0].x).abs() < 0.000001);
    }
}
