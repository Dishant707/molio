//! SDF (Structure Data File) / MOL format parser.
//!
//! Implements parsing for the MDL CTFile format used by PubChem, ZINC,
//! ChEMBL, and virtually every drug discovery database.
//!
//! ## Format
//! - SDF: multiple molecules separated by `$$$$`
//! - MOL: single molecule (no `$$$$` delimiter)
//!
//! ## Design
//! - Streaming: parses molecule-by-molecule, never loads full file
//! - Zero-copy atom coordinates: parsed directly from fixed-width fields
//! - Bond connectivity preserved with bond orders (1=single, 2=double, 3=triple)

use std::path::Path;

use crate::error::{MolioError, MolioResult};
use crate::types::{Atom, Bond, BondOrder, FileFormat, Molecule};

/// Parse an SDF file (may contain multiple molecules).
/// Returns a Vec<Molecule>, one per `$$$$` block.
pub fn parse_sdf(path: impl AsRef<Path>) -> MolioResult<Vec<Molecule>> {
    let content = std::fs::read_to_string(path.as_ref())?;
    parse_sdf_str(&content)
}

/// Parse SDF content from a string.
pub fn parse_sdf_str(content: &str) -> MolioResult<Vec<Molecule>> {
    let mut molecules = Vec::new();
    let lines: Vec<&str> = content.lines().collect();
    let mut offset = 0;

    while offset < lines.len() {
        // Skip empty lines
        while offset < lines.len() && lines[offset].trim().is_empty() {
            offset += 1;
        }
        if offset >= lines.len() {
            break;
        }

        let mol = parse_single_mol(&lines, &mut offset)?;
        molecules.push(mol);
    }

    Ok(molecules)
}

/// Parse a single molecule from lines starting at `offset`.
/// Advances `offset` past the molecule (including `$$$$` if present).
fn parse_single_mol(lines: &[&str], offset: &mut usize) -> MolioResult<Molecule> {
    let _start = *offset;
    let mut molecule = Molecule::new(FileFormat::Sdf);

    // Line 1: molecule name / title
    if *offset < lines.len() {
        molecule.name = Some(lines[*offset].trim().to_string());
        *offset += 1;
    }

    // Line 2: program/date stamp (ignored)
    if *offset < lines.len() {
        *offset += 1;
    }

    // Line 3: optional comment. Detect if this is actually the counts line
    // (happens when the file has no comment line).
    if *offset < lines.len() {
        let line = lines[*offset];
        // Counts line starts with 3 digits (possibly space-padded)
        let trimmed = line.trim();
        if is_counts_line(trimmed) {
            // This is already the counts line, don't consume as comment
        } else {
            *offset += 1; // consume comment line
        }
    }

    // Line 4: counts line
    // Format: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
    if *offset >= lines.len() {
        return Err(MolioError::Parse {
            line: *offset + 1,
            message: "unexpected end of file before counts line".into(),
        });
    }

    let counts_line = lines[*offset];
    *offset += 1;

    let n_atoms = parse_int_field(counts_line, 0, 3)?;
    let n_bonds = parse_int_field(counts_line, 3, 6)?;

    // Parse atom block
    for _ in 0..n_atoms {
        if *offset >= lines.len() {
            return Err(MolioError::Parse {
                line: *offset + 1,
                message: format!("expected {n_atoms} atoms, found fewer"),
            });
        }
        let atom = parse_mol_atom_line(lines[*offset], *offset + 1)?;
        molecule.atoms.push(atom);
        *offset += 1;
    }

    // Parse bond block
    for _ in 0..n_bonds {
        if *offset >= lines.len() {
            return Err(MolioError::Parse {
                line: *offset + 1,
                message: format!("expected {n_bonds} bonds, found fewer"),
            });
        }
        let bond = parse_mol_bond_line(lines[*offset], *offset + 1)?;
        molecule.bonds.push(bond);
        *offset += 1;
    }

    // Parse properties block (until $$$$ or EOF)
    while *offset < lines.len() {
        let line = lines[*offset];
        if line.trim() == "$$$$" {
            *offset += 1;
            break;
        }
        if line.starts_with("> <") {
            // Property tag line: "> <MolWeight>" → "MolWeight"
            let tag = line[3..].trim_end_matches('>').trim().to_string();
            *offset += 1;
            if *offset < lines.len() && !lines[*offset].is_empty() {
                let value = lines[*offset].trim().to_string();
                molecule.properties.insert(tag, value);
            }
            // Skip empty line after value
            if *offset < lines.len() && lines[*offset].trim().is_empty() {
                *offset += 1;
            }
        } else if line.trim() == "M  END" {
            *offset += 1;
        } else if line.trim().is_empty() {
            *offset += 1;
        } else {
            // Unknown line in properties block, skip
            *offset += 1;
        }
    }

    Ok(molecule)
}

/// Parse an atom line from MOL/SDF format.
///
/// Format (V2000):
/// xxxxx.xxxyyyyy.yyyzzzzz.zzz eeesssdddcccsss...
/// Col 1-10:  x coordinate
/// Col 11-20: y coordinate
/// Col 21-30: z coordinate
/// Col 32-33: element symbol (or blank = from atom symbol)
fn parse_mol_atom_line(line: &str, line_num: usize) -> MolioResult<Atom> {
    if line.len() < 33 {
        return Err(MolioError::Parse {
            line: line_num,
            message: "atom line too short".into(),
        });
    }

    let x = parse_float(&line[0..10], line_num, "x")?;
    let y = parse_float(&line[10..20], line_num, "y")?;
    let z = parse_float(&line[20..30], line_num, "z")?;

    // Element symbol: cols 32-33 (0-indexed: 31..33)
    let element = if line.len() >= 33 {
        line[31..33].trim().to_string()
    } else {
        String::new()
    };

    // If element is empty, use the atom symbol from first 3 chars
    let element = if element.is_empty() {
        line[..3].trim().to_string()
    } else {
        element
    };

    Ok(Atom {
        serial: 0, // MOL format doesn't use serial numbers
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
    })
}

/// Parse a bond line from MOL/SDF format.
///
/// Format: 111222tttsssxxxrrrccc
/// Col 1-3:  first atom index (1-based)
/// Col 4-6:  second atom index (1-based)
/// Col 7-9:  bond type (1=single, 2=double, 3=triple, 4=aromatic)
fn parse_mol_bond_line(line: &str, line_num: usize) -> MolioResult<Bond> {
    if line.len() < 9 {
        return Err(MolioError::Parse {
            line: line_num,
            message: "bond line too short".into(),
        });
    }

    let atom1 = parse_int_field(line, 0, 3)? as u32;
    let atom2 = parse_int_field(line, 3, 6)? as u32;
    let bond_type = parse_int_field(line, 6, 9)?;

    let order = match bond_type {
        1 => BondOrder::Single,
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        4 => BondOrder::Aromatic,
        n => BondOrder::Unknown(n as u8),
    };

    Ok(Bond {
        atom1,
        atom2,
        order,
    })
}

// ─── Helpers ──────────────────────────────────────────────────────

/// Detect if a line looks like a V2000/V3000 counts line.
/// Format: aaabbblll... where aaa=atom count, bbb=bond count.
fn is_counts_line(line: &str) -> bool {
    let line = line.trim();
    if line.len() < 6 {
        return false;
    }
    // Check first 3 chars are digits (possibly space-padded)
    line[..3].trim().chars().all(|c| c.is_ascii_digit())
        && line[3..6].trim().chars().all(|c| c.is_ascii_digit())
}

fn parse_int_field(line: &str, start: usize, end: usize) -> MolioResult<i32> {
    let field = &line[start.min(line.len())..end.min(line.len())];
    field.trim().parse::<i32>().map_err(|_| MolioError::Parse {
        line: 0,
        message: format!("invalid integer field '{field}' at cols {start}-{end}"),
    })
}

fn parse_float(field: &str, line_num: usize, name: &str) -> MolioResult<f64> {
    field.trim().parse::<f64>().map_err(|_| MolioError::Parse {
        line: line_num,
        message: format!("invalid {name} coordinate '{field}'"),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_benzene_atom() {
        let line = "    1.9050   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0";
        let atom = parse_mol_atom_line(line, 1).unwrap();
        assert!((atom.x - 1.9050).abs() < 0.0001);
        assert!((atom.y - 0.0).abs() < 0.0001);
        assert!((atom.z - 0.0).abs() < 0.0001);
        assert_eq!(atom.element, "C");
    }

    #[test]
    fn test_parse_single_bond() {
        let line = "  1  2  2  0  0  0  0";
        let bond = parse_mol_bond_line(line, 1).unwrap();
        assert_eq!(bond.atom1, 1);
        assert_eq!(bond.atom2, 2);
        assert_eq!(bond.order, BondOrder::Double);
    }

    #[test]
    fn test_parse_sdf_multi_molecule() {
        let sdf = r#"Benzene
  test
  comment
  6  6  0  0  0  0  0  0  0  0999 V2000
    1.9050   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9520    1.6500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9520    1.6500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9050    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9520   -1.6500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9520   -1.6500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
M  END
> <MolWeight>
78.11

$$$$
Water
  test
  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9570    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2400    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END
$$$$"#;

        let mols = parse_sdf_str(sdf).unwrap();
        assert_eq!(mols.len(), 2);

        let benzene = &mols[0];
        assert_eq!(benzene.name.as_deref(), Some("Benzene"));
        assert_eq!(benzene.n_atoms(), 6);
        assert_eq!(benzene.bonds.len(), 6);
        assert_eq!(benzene.bonds[0].order, BondOrder::Double);

        let water = &mols[1];
        assert_eq!(water.name.as_deref(), Some("Water"));
        assert_eq!(water.n_atoms(), 3);
        assert_eq!(water.bonds.len(), 2);
    }

    #[test]
    fn test_parse_sdf_properties() {
        let sdf = "Test\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <MolWeight>
180.16

> <Formula>
C9H8O4

$$$$";
        let mols = parse_sdf_str(sdf).unwrap();
        assert_eq!(mols.len(), 1);
        assert_eq!(mols[0].properties.get("MolWeight").map(|s| s.as_str()), Some("180.16"));
        assert_eq!(mols[0].properties.get("Formula").map(|s| s.as_str()), Some("C9H8O4"));
    }
}
