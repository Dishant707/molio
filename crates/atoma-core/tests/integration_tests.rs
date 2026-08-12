//! Integration tests for molio-core.
//!
//! These tests verify correctness by cross-validating against
//! known reference data and testing roundtrip fidelity.

use atoma_core::parser::pdb::parse_pdb_str;
use atoma_core::parser::sdf::parse_sdf_str;

const PDB_1CRN: &str = include_str!("../../../test_data/pdb/1crn.pdb");

#[test]
fn test_parse_1crn_atom_count() {
    let mol = parse_pdb_str(PDB_1CRN, "1crn.pdb").expect("should parse 1CRN");
    // 1CRN residues 1-10 (first 10 residues of crambin)
    assert_eq!(mol.n_atoms(), 71, "1CRN should have 71 atoms");
}

#[test]
fn test_parse_1crn_chains() {
    let mol = parse_pdb_str(PDB_1CRN, "1crn.pdb").expect("should parse 1CRN");
    assert_eq!(mol.chains.len(), 1, "1CRN should have 1 chain");
    assert_eq!(mol.chains[0].id, 'A');
}

#[test]
fn test_parse_1crn_residues() {
    let mol = parse_pdb_str(PDB_1CRN, "1crn.pdb").expect("should parse 1CRN");
    assert_eq!(mol.n_residues(), 10, "1CRN should have 10 residues");

    let expected_residues = ["THR", "THR", "CYS", "CYS", "PRO", "SER", "ILE", "VAL", "ARG", "SER"];
    for (i, res) in mol.chains[0].residues.iter().enumerate() {
        assert_eq!(res.name, expected_residues[i], "Residue {} should be {}", i + 1, expected_residues[i]);
    }
}

#[test]
fn test_parse_1crn_first_atom() {
    let mol = parse_pdb_str(PDB_1CRN, "1crn.pdb").expect("should parse 1CRN");
    let first = &mol.atoms[0];

    assert_eq!(first.serial, 1);
    assert_eq!(first.name, "N");
    assert_eq!(first.res_name, "THR");
    assert_eq!(first.chain_id, 'A');
    assert_eq!(first.res_seq, 1);
    assert!((first.x - 9.670).abs() < 0.001);
    assert!((first.y - 10.289).abs() < 0.001);
    assert!((first.z - 11.135).abs() < 0.001);
    assert_eq!(first.element, "N");
}

#[test]
fn test_parse_1crn_last_atom() {
    let mol = parse_pdb_str(PDB_1CRN, "1crn.pdb").expect("should parse 1CRN");
    let last = mol.atoms.last().unwrap();

    assert_eq!(last.serial, 71);
    assert_eq!(last.name, "OG");
    assert_eq!(last.res_name, "SER");
    assert_eq!(last.res_seq, 10);
    assert_eq!(last.element, "O");
}

#[test]
fn test_roundtrip_coordinate_precision() {
    // Parse 1CRN twice and verify identical results
    let mol1 = parse_pdb_str(PDB_1CRN, "1crn.pdb").expect("first parse");
    let mol2 = parse_pdb_str(PDB_1CRN, "1crn.pdb").expect("second parse");

    assert_eq!(mol1.n_atoms(), mol2.n_atoms());
    for (a1, a2) in mol1.atoms.iter().zip(mol2.atoms.iter()) {
        assert_eq!(a1.serial, a2.serial);
        assert_eq!(a1.x, a2.x);
        assert_eq!(a1.y, a2.y);
        assert_eq!(a1.z, a2.z);
        assert_eq!(a1.element, a2.element);
    }
}

#[test]
fn test_molecular_weight() {
    let mol = parse_pdb_str(PDB_1CRN, "1crn.pdb").expect("should parse 1CRN");
    let mw = mol.molecular_weight();
    // Should be ~900 Da for 10 small residues
    assert!(mw > 800.0 && mw < 1200.0,
        "Molecular weight {} should be in range 800-1200", mw);
}

#[test]
fn test_empty_input() {
    let result = parse_pdb_str("", "empty.pdb");
    assert!(result.is_ok());
    let mol = result.unwrap();
    assert_eq!(mol.n_atoms(), 0);
}

#[test]
fn test_header_parsing() {
    let mol = parse_pdb_str(PDB_1CRN, "1crn.pdb").expect("should parse");
    assert!(!mol.header.is_empty());
    assert!(mol.header[0].contains("HEADER"));
}

#[test]
fn test_alternate_conformation() {
    // Atom with alternate location marker
    let line = "ATOM     10  CA AALA A   2      12.345   0.000   0.000  0.50 20.00           C  ";
    let mol = parse_pdb_str(line, "test.pdb").unwrap();
    assert_eq!(mol.atoms[0].alt_loc, 'A');
    assert_eq!(mol.atoms[0].occupancy, 0.50);
}

#[test]
fn test_negative_coordinates() {
    let pdb = "\
ATOM      1  N   ALA A   1     -12.345 -67.890  -0.500  1.00  0.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
";
    let mol = parse_pdb_str(pdb, "test.pdb").unwrap();
    assert_eq!(mol.atoms[0].x, -12.345);
    assert_eq!(mol.atoms[0].y, -67.890);
    assert_eq!(mol.atoms[0].z, -0.500);
}

#[test]
fn test_hetatm_parsing() {
    let pdb = "\
HETATM 1001  O   HOH A 101      15.000  20.000  25.000  1.00 30.00           O
HETATM 1002  ZN  ZN  A 102      10.000  10.000  10.000  1.00 15.00          ZN
";
    let mol = parse_pdb_str(pdb, "test.pdb").unwrap();
    assert_eq!(mol.n_atoms(), 2);
    assert_eq!(mol.atoms[0].res_name, "HOH");
    assert_eq!(mol.atoms[1].element, "ZN");
}

#[test]
fn test_bounding_box() {
    let mol = parse_pdb_str(PDB_1CRN, "1crn.pdb").unwrap();
    let (min, max) = mol.bounding_box().expect("should have bounding box");
    // 1CRN fragment (10 residues) is ~18Å across in X/Y
    assert!(max[0] - min[0] > 5.0, "X range too small: {}", max[0] - min[0]);
    assert!(max[1] - min[1] > 5.0, "Y range too small: {}", max[1] - min[1]);
    assert!(max[2] - min[2] > 2.0, "Z range too small: {}", max[2] - min[2]);
}

// ─── SDF Integration Tests ────────────────────────────────────────

const SDF_MOLECULES: &str = include_str!("../../../test_data/sdf/molecules.sdf");

#[test]
fn test_parse_sdf_three_molecules() {
    let mols = parse_sdf_str(SDF_MOLECULES).expect("should parse SDF");
    assert_eq!(mols.len(), 3, "should have 3 molecules");
}

#[test]
fn test_parse_sdf_aspirin() {
    let mols = parse_sdf_str(SDF_MOLECULES).expect("should parse SDF");
    let aspirin = &mols[0];
    assert_eq!(aspirin.name.as_deref(), Some("Aspirin"));
    assert_eq!(aspirin.n_atoms(), 21);
    assert_eq!(aspirin.bonds.len(), 21);
    assert_eq!(aspirin.properties.get("Formula").map(|s| s.as_str()), Some("C9H8O4"));
    assert_eq!(aspirin.properties.get("MolWeight").map(|s| s.as_str()), Some("180.16"));
}

#[test]
fn test_parse_sdf_caffeine() {
    let mols = parse_sdf_str(SDF_MOLECULES).expect("should parse SDF");
    let caffeine = &mols[1];
    assert_eq!(caffeine.name.as_deref(), Some("Caffeine"));
    assert_eq!(caffeine.n_atoms(), 14);
    assert_eq!(caffeine.bonds.len(), 15);
}

#[test]
fn test_parse_sdf_benzene() {
    let mols = parse_sdf_str(SDF_MOLECULES).expect("should parse SDF");
    let benzene = &mols[2];
    assert_eq!(benzene.name.as_deref(), Some("Benzene"));
    assert_eq!(benzene.n_atoms(), 6);
    assert_eq!(benzene.bonds.len(), 6);
}

#[test]
fn test_parse_sdf_bond_orders() {
    let mols = parse_sdf_str(SDF_MOLECULES).expect("should parse SDF");
    let benzene = &mols[2];
    // Benzene: 3 double bonds, 3 single bonds (Kekule form)
    let doubles = benzene.bonds.iter().filter(|b| b.order == atoma_core::BondOrder::Double).count();
    let singles = benzene.bonds.iter().filter(|b| b.order == atoma_core::BondOrder::Single).count();
    assert_eq!(doubles, 3);
    assert_eq!(singles, 3);
}

#[test]
fn test_parse_sdf_coordinates() {
    let mols = parse_sdf_str(SDF_MOLECULES).expect("should parse SDF");
    let benzene = &mols[2];
    // First carbon should be at (1.905, 0, 0)
    let c1 = &benzene.atoms[0];
    assert!((c1.x - 1.9050).abs() < 0.001);
    assert!((c1.y - 0.0).abs() < 0.001);
    assert_eq!(c1.element, "C");
}

#[test]
fn test_parse_sdf_roundtrip() {
    let mols1 = parse_sdf_str(SDF_MOLECULES).expect("first parse");
    let mols2 = parse_sdf_str(SDF_MOLECULES).expect("second parse");
    assert_eq!(mols1.len(), mols2.len());
    for (m1, m2) in mols1.iter().zip(mols2.iter()) {
        assert_eq!(m1.n_atoms(), m2.n_atoms());
        assert_eq!(m1.bonds.len(), m2.bonds.len());
        for (a1, a2) in m1.atoms.iter().zip(m2.atoms.iter()) {
            assert_eq!(a1.x, a2.x);
            assert_eq!(a1.y, a2.y);
            assert_eq!(a1.element, a2.element);
        }
    }
}
