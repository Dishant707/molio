//! XYZ format parser.
//!
//! XYZ is the simplest molecular format — an atom-count header line,
//! a comment line, then one line per atom: element x y z.
//!
//! Used by most computational chemistry software (Gaussian, ORCA, VASP).

use std::path::Path;

use crate::error::{MolioError, MolioResult};
use crate::types::{Atom, FileFormat, Molecule};

pub fn parse_xyz(path: impl AsRef<Path>) -> MolioResult<Molecule> {
    let content = std::fs::read_to_string(path.as_ref())?;
    parse_xyz_str(&content)
}

pub fn parse_xyz_str(content: &str) -> MolioResult<Molecule> {
    let mut molecule = Molecule::new(FileFormat::Xyz);
    let mut lines = content.lines().enumerate();

    // Line 1: atom count
    let (_, first_line) = lines.next().ok_or_else(|| MolioError::Parse {
        line: 1,
        message: "empty XYZ file".into(),
    })?;
    let n_atoms: usize = first_line.trim().parse().map_err(|_| MolioError::Parse {
        line: 1,
        message: format!("invalid atom count: '{}'", first_line.trim()),
    })?;

    // Line 2: comment/title (optional)
    if let Some((_, comment)) = lines.next() {
        molecule.name = Some(comment.trim().to_string());
    }

    // Parse atom lines
    let mut serial = 0u32;
    for (line_num, line) in lines {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            continue; // Skip malformed lines
        }

        let element = parts[0].to_string();
        let x: f64 = parts[1].parse().map_err(|_| MolioError::Parse {
            line: line_num + 1,
            message: format!("invalid x coordinate: '{}'", parts[1]),
        })?;
        let y: f64 = parts[2].parse().map_err(|_| MolioError::Parse {
            line: line_num + 1,
            message: format!("invalid y coordinate: '{}'", parts[2]),
        })?;
        let z: f64 = parts[3].parse().map_err(|_| MolioError::Parse {
            line: line_num + 1,
            message: format!("invalid z coordinate: '{}'", parts[3]),
        })?;

        serial += 1;
        molecule.atoms.push(Atom {
            serial,
            name: element.clone(),
            alt_loc: ' ',
            res_name: String::new(),
            chain_id: ' ',
            res_seq: 0,
            i_code: ' ',
            x,
            y,
            z,
            occupancy: 1.0,
            temp_factor: 0.0,
            element,
            charge: None,
        });
    }

    // Validate atom count matches
    if molecule.n_atoms() != n_atoms {
        return Err(MolioError::Parse {
            line: 1,
            message: format!(
                "atom count mismatch: header says {}, found {}",
                n_atoms,
                molecule.n_atoms()
            ),
        });
    }

    Ok(molecule)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_water() {
        let xyz = "3\nWater\nO  0.000  0.000  0.000\nH  0.957  0.000  0.000\nH -0.240  0.927  0.000\n";
        let mol = parse_xyz_str(xyz).unwrap();
        assert_eq!(mol.n_atoms(), 3);
        assert_eq!(mol.name.as_deref(), Some("Water"));
        assert_eq!(mol.atoms[0].element, "O");
        assert_eq!(mol.atoms[1].element, "H");
    }

    #[test]
    fn test_parse_no_comment() {
        let xyz = "1\n\nC  0.0  0.0  0.0\n";
        let mol = parse_xyz_str(xyz).unwrap();
        assert_eq!(mol.n_atoms(), 1);
    }

    #[test]
    fn test_count_mismatch() {
        let xyz = "5\nTest\nC  0.0  0.0  0.0\n";
        assert!(parse_xyz_str(xyz).is_err());
    }

    #[test]
    fn test_roundtrip() {
        let xyz = "2\nEthane\nC  0.0  0.0  0.0\nC  1.5  0.0  0.0\n";
        let m1 = parse_xyz_str(xyz).unwrap();
        let m2 = parse_xyz_str(xyz).unwrap();
        assert_eq!(m1.n_atoms(), m2.n_atoms());
        for (a1, a2) in m1.atoms.iter().zip(m2.atoms.iter()) {
            assert_eq!(a1.x, a2.x);
            assert_eq!(a1.element, a2.element);
        }
    }
}
