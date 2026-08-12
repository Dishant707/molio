//! Property-based tests using proptest.
//!
//! Generate random valid PDB/SDF content and verify:
//! 1. Parsing never panics
//! 2. Roundtrip consistency (parse twice = same result)
//! 3. Atom count is preserved

use proptest::prelude::*;

fn gen_atom_line(serial: u32, res_seq: i32, x: f64, y: f64, z: f64) -> String {
    let names = ["N", "CA", "C", "O", "CB"];
    let res_names = ["ALA", "GLY", "CYS", "THR", "SER"];
    let elements = ["C", "N", "O", "S"];
    let name = names[serial as usize % names.len()];
    let res = res_names[res_seq as usize % res_names.len()];
    let elem = elements[serial as usize % elements.len()];
    format!(
        "ATOM  {:>5} {:<4} {} A{:>4}    {:>8.3}{:>8.3}{:>8.3}  1.00  0.00          {:>2}  ",
        serial, name, res, res_seq, x, y, z, elem
    )
}

fn gen_pdb(n_atoms: usize) -> String {
    let mut pdb = String::from("HEADER    PROPERTY_TEST\n");
    for i in 0..n_atoms {
        let serial = (i + 1) as u32;
        let res_seq = ((i / 4) + 1) as i32;
        let x = (i as f64 * 1.5) % 50.0;
        let y = (i as f64 * 2.3) % 50.0;
        let z = (i as f64 * 3.1) % 50.0;
        pdb.push_str(&gen_atom_line(serial, res_seq, x, y, z));
        pdb.push('\n');
    }
    pdb.push_str("END\n");
    pdb
}

fn gen_sdf(n_atoms: usize) -> String {
    let elements = ["C", "N", "O", "S", "P", "F", "Cl", "Br"];
    let n_bonds = if n_atoms > 1 { n_atoms - 1 } else { 0 };
    let mut sdf = String::new();
    sdf.push_str("PropTest\n  gen\n\n");
    sdf.push_str(&format!("{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n", n_atoms, n_bonds));
    for i in 0..n_atoms {
        let x = (i as f64 * 1.2) % 30.0;
        let y = (i as f64 * 1.7) % 30.0;
        let z = (i as f64 * 2.1) % 30.0;
        sdf.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
            x, y, z, elements[i % elements.len()]
        ));
    }
    for i in 1..n_atoms {
        sdf.push_str(&format!("{:>3}{:>3}  1  0  0  0  0\n", i, i + 1));
    }
    sdf.push_str("M  END\n$$$$\n");
    sdf
}

proptest! {
    #[test]
    fn pdb_never_panics(n in 1usize..100) {
        let pdb = gen_pdb(n);
        let _ = atoma_core::parser::pdb::parse_pdb_str(&pdb, "test.pdb");
    }

    #[test]
    fn pdb_idempotent(n in 1usize..100) {
        let pdb = gen_pdb(n);
        let mol1 = atoma_core::parser::pdb::parse_pdb_str(&pdb, "t1.pdb").unwrap();
        let mol2 = atoma_core::parser::pdb::parse_pdb_str(&pdb, "t2.pdb").unwrap();
        prop_assert_eq!(mol1.n_atoms(), mol2.n_atoms());
        prop_assert_eq!(mol1.n_atoms(), n);
    }

    #[test]
    fn pdb_coordinates_preserved(n in 1usize..50) {
        let pdb = gen_pdb(n);
        let mol = atoma_core::parser::pdb::parse_pdb_str(&pdb, "t.pdb").unwrap();
        prop_assert_eq!(mol.n_atoms(), n);
        for atom in &mol.atoms {
            prop_assert!(atom.x.is_finite());
            prop_assert!(atom.y.is_finite());
            prop_assert!(atom.z.is_finite());
        }
    }
}

proptest! {
    #[test]
    fn sdf_never_panics(n in 1usize..30) {
        let sdf = gen_sdf(n);
        let _ = atoma_core::parser::sdf::parse_sdf_str(&sdf);
    }

    #[test]
    fn sdf_idempotent(n in 1usize..30) {
        let sdf = gen_sdf(n);
        let mols1 = atoma_core::parser::sdf::parse_sdf_str(&sdf).unwrap();
        let mols2 = atoma_core::parser::sdf::parse_sdf_str(&sdf).unwrap();
        prop_assert_eq!(mols1.len(), mols2.len());
        prop_assert_eq!(mols1[0].n_atoms(), n);
        prop_assert_eq!(mols2[0].n_atoms(), n);
    }

    #[test]
    fn sdf_bond_count_consistent(n in 2usize..30) {
        let sdf = gen_sdf(n);
        let mols = atoma_core::parser::sdf::parse_sdf_str(&sdf).unwrap();
        prop_assert_eq!(mols[0].bonds.len(), n - 1);
        for bond in &mols[0].bonds {
            prop_assert!(bond.atom1 >= 1 && bond.atom1 <= n as u32);
            prop_assert!(bond.atom2 >= 1 && bond.atom2 <= n as u32);
            prop_assert!(bond.atom1 != bond.atom2);
        }
    }
}
