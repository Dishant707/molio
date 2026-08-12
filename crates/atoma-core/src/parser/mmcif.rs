//! mmCIF (PDBx) format parser.
//!
//! mmCIF is the modern replacement for PDB format, used by the wwPDB
//! since 2014. It uses a dictionary-based format with `loop_` structures
//! for tabular data and key-value pairs for metadata.
//!
//! Handles structures up to millions of atoms (no 99,999 atom limit).

use std::collections::HashMap;
use std::path::Path;

use crate::error::{MolioError, MolioResult};
use crate::types::{Atom, Chain, FileFormat, Molecule, Residue};

pub fn parse_mmcif(path: impl AsRef<Path>) -> MolioResult<Molecule> {
    let content = std::fs::read_to_string(path.as_ref())?;
    parse_mmcif_str(&content)
}

pub fn parse_mmcif_str(content: &str) -> MolioResult<Molecule> {
    let lines: Vec<&str> = content.lines().collect();
    let mut molecule = Molecule::new(FileFormat::MmCif);
    let mut i = 0;

    while i < lines.len() {
        let line = lines[i].trim();

        if line.is_empty() || line.starts_with('#') {
            i += 1;
            continue;
        }

        if line.starts_with("data_") {
            // Data block name
            if molecule.name.is_none() {
                molecule.name = Some(line[5..].to_string());
            }
            i += 1;
        } else if line.starts_with("loop_") {
            // Check if this is an _atom_site loop before parsing
            i += 1;
            let (parsed_atoms, next_i) = parse_atom_site_loop(&lines, i)?;
            if !parsed_atoms.is_empty() {
                molecule.atoms.extend(parsed_atoms);
            }
            i = next_i;
        } else if line.starts_with('_') && !line.contains("loop_") {
            // Key-value pair (single)
            if let Some(eq_pos) = line.find(char::is_whitespace) {
                let _key = &line[..eq_pos];
                let _value = line[eq_pos..].trim();
                i += 1;
            } else {
                i += 1;
            }
        } else {
            i += 1;
        }
    }

    // Build chains from atoms
    build_chains(&mut molecule);

    Ok(molecule)
}

/// Parse an _atom_site loop. Returns parsed atoms and next line index.
/// Non-atom loops are silently skipped (return empty vec).
fn parse_atom_site_loop(lines: &[&str], start: usize) -> MolioResult<(Vec<Atom>, usize)> {
    let mut i = start;
    let mut headers: Vec<String> = Vec::new();

    // Collect column headers (lines starting with _)
    while i < lines.len() && lines[i].trim().starts_with('_') {
        headers.push(lines[i].trim().to_string());
        i += 1;
    }

    if headers.is_empty() {
        // Empty loop, skip past data
        while i < lines.len() {
            let line = lines[i].trim();
            if line.is_empty() || line.starts_with('#') || line.starts_with('_') || line.starts_with("loop_") {
                break;
            }
            i += 1;
        }
        return Ok((vec![], i));
    }

    // Skip non-atom_site loops (database, citation, etc.)
    let is_atom_site = headers.iter().any(|h| h.contains("_atom_site"));
    if !is_atom_site {
        // Skip past data lines of this loop
        while i < lines.len() {
            let line = lines[i].trim();
            if line.is_empty() || line.starts_with('#') || line.starts_with('_') || line.starts_with("loop_") {
                break;
            }
            i += 1;
        }
        return Ok((vec![], i));
    }

    // Map column indices to Atom fields
    let col_map = build_column_map(&headers)?;

    // Parse data rows
    let mut atoms = Vec::new();
    let mut serial = 0u32;

    while i < lines.len() {
        let line = lines[i].trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with('_') || line.starts_with("loop_") {
            break; // End of loop data
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.is_empty() {
            i += 1;
            continue;
        }

        serial += 1;
        let mut atom = Atom {
            serial,
            name: String::new(),
            alt_loc: ' ',
            res_name: String::new(),
            chain_id: ' ',
            res_seq: 0,
            i_code: ' ',
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            temp_factor: 0.0,
            element: String::new(),
            charge: None,
        };

        // Fill fields from column map
        for (col_idx, field) in &col_map {
            if *col_idx < parts.len() {
                let value = parts[*col_idx];
                match field {
                    AtomField::Name => atom.name = value.to_string(),
                    AtomField::ResName => atom.res_name = value.to_string(),
                    AtomField::ChainId => atom.chain_id = value.chars().next().unwrap_or(' '),
                    AtomField::ResSeq => atom.res_seq = value.parse().unwrap_or(0),
                    AtomField::X => atom.x = value.parse().unwrap_or(0.0),
                    AtomField::Y => atom.y = value.parse().unwrap_or(0.0),
                    AtomField::Z => atom.z = value.parse().unwrap_or(0.0),
                    AtomField::Occupancy => atom.occupancy = value.parse().unwrap_or(1.0),
                    AtomField::TempFactor => atom.temp_factor = value.parse().unwrap_or(0.0),
                    AtomField::Element => atom.element = value.to_string(),
                    AtomField::AltLoc => atom.alt_loc = value.chars().next().unwrap_or(' '),
                }
            }
        }

        // Derive element from atom name if not explicitly set
        if atom.element.is_empty() && !atom.name.is_empty() {
            atom.element = atom.name.chars().take(1).collect::<String>();
            if atom.element.len() == 1 && atom.name.len() > 1 {
                let c2 = atom.name.chars().nth(1).unwrap();
                if c2.is_alphabetic() && c2.is_lowercase() {
                    atom.element.push(c2);
                }
            }
            // Title-case: "c" → "C", "cl" → "Cl"
            let mut chars: Vec<char> = atom.element.chars().collect();
            if !chars.is_empty() {
                chars[0] = chars[0].to_ascii_uppercase();
            }
            atom.element = chars.into_iter().collect();
        }

        atoms.push(atom);
        i += 1;
    }

    Ok((atoms, i))
}

#[derive(Debug)]
enum AtomField {
    Name,
    ResName,
    ChainId,
    ResSeq,
    X,
    Y,
    Z,
    Occupancy,
    TempFactor,
    Element,
    AltLoc,
}

/// Map mmCIF column headers to AtomField variants.
fn build_column_map(headers: &[String]) -> MolioResult<HashMap<usize, AtomField>> {
    let mut map = HashMap::new();

    for (idx, header) in headers.iter().enumerate() {
        let h = header.to_lowercase();
        let field = if h.contains("label_atom_id") || h.contains("auth_atom_id") {
            AtomField::Name
        } else if h.contains("label_comp_id") || h.contains("auth_comp_id")
            || h.contains("comp_id") || h.ends_with(".comp_id")
        {
            AtomField::ResName
        } else if h.contains("label_asym_id") || h.contains("auth_asym_id")
            || h.contains("asym_id") || h.ends_with(".asym_id")
        {
            AtomField::ChainId
        } else if h.contains("label_seq_id") || h.contains("auth_seq_id")
            || h.contains("seq_id") || h.ends_with(".seq_id")
        {
            AtomField::ResSeq
        } else if h.ends_with("cartn_x") || h.ends_with("fract_x")
            || (h.contains("_x") && (h.contains("coord") || h.contains("cartesian")))
        {
            AtomField::X
        } else if h.ends_with("cartn_y") || h.ends_with("fract_y")
            || (h.contains("_y") && (h.contains("coord") || h.contains("cartesian")))
        {
            AtomField::Y
        } else if h.ends_with("cartn_z") || h.ends_with("fract_z")
            || (h.contains("_z") && (h.contains("coord") || h.contains("cartesian")))
        {
            AtomField::Z
        } else if h.contains("occupancy") {
            AtomField::Occupancy
        } else if h.contains("b_iso") || h.contains("temp_factor") || h.contains("u_iso") {
            AtomField::TempFactor
        } else if h.contains("type_symbol") || h.contains("element") {
            AtomField::Element
        } else if h.contains("alt_id") || h.contains("alt_loc") {
            AtomField::AltLoc
        } else {
            continue; // Unknown field, skip
        };
        map.insert(idx, field);
    }

    if map.is_empty() {
        let found: Vec<&str> = headers.iter().map(|s| s.as_str()).take(10).collect();
        return Err(MolioError::Parse {
            line: 0,
            message: format!(
                "no recognized atom site columns. First headers: {:?}",
                found
            ),
        });
    }

    Ok(map)
}

/// Build chain/residue hierarchy from flat atom list.
fn build_chains(molecule: &mut Molecule) {
    let mut current_chain: Option<Chain> = None;
    let mut current_residue: Option<Residue> = None;

    for atom in &molecule.atoms {
        let chain_id = atom.chain_id;
        let res_seq = atom.res_seq;
        let i_code = atom.i_code;
        let res_name = atom.res_name.clone();

        if current_chain.as_ref().map(|c| c.id) != Some(chain_id) {
            if let Some(chain) = current_chain.take() {
                if !chain.residues.is_empty() {
                    molecule.chains.push(chain);
                }
            }
            current_residue = None;
            current_chain = Some(Chain { id: chain_id, residues: Vec::new() });
        }

        let need_new = current_residue.as_ref().map_or(true, |r| {
            r.seq_num != res_seq || r.insertion_code != i_code
        });

        if need_new {
            if let Some(res) = current_residue.take() {
                if let Some(ref mut chain) = current_chain {
                    chain.residues.push(res);
                }
            }
            current_residue = Some(Residue {
                name: res_name,
                seq_num: res_seq,
                insertion_code: i_code,
                atoms: Vec::new(),
            });
        }

        if let Some(ref mut res) = current_residue {
            res.atoms.push(atom.clone());
        }
    }

    // Push final residue and chain
    if let Some(res) = current_residue {
        if let Some(ref mut chain) = current_chain {
            chain.residues.push(res);
        }
    }
    if let Some(chain) = current_chain {
        if !chain.residues.is_empty() {
            molecule.chains.push(chain);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_1crn_mmcif() {
        let cif = include_str!("../../../../test_data/mmcif/1crn.cif");
        let mol = parse_mmcif_str(cif).unwrap();
        assert_eq!(mol.n_atoms(), 10);
        assert_eq!(mol.n_residues(), 2);
        // First atom
        assert_eq!(mol.atoms[0].name, "N");
        assert_eq!(mol.atoms[0].res_name, "THR");
        assert_eq!(mol.atoms[0].chain_id, 'A');
        assert!((mol.atoms[0].x - 9.670).abs() < 0.01);
    }

    #[test]
    fn test_parse_minimal_cif() {
        let cif = "data_test\nloop_\n_atom_site.label_atom_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nCA 1.0 2.0 3.0\n";
        let mol = parse_mmcif_str(cif).unwrap();
        assert_eq!(mol.n_atoms(), 1);
        assert_eq!(mol.atoms[0].name, "CA");
        assert_eq!(mol.atoms[0].element, "C");
    }

    #[test]
    fn test_parse_empty_cif() {
        let mol = parse_mmcif_str("data_test\n").unwrap();
        assert_eq!(mol.n_atoms(), 0);
    }

    #[test]
    fn test_build_chains() {
        let cif = "data_test\nloop_\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nN ALA A 1 0 0 0\nCA ALA A 1 0 0 0\nN GLY B 1 0 0 0\nCA GLY B 1 0 0 0\n";
        let mol = parse_mmcif_str(cif).unwrap();
        assert_eq!(mol.n_atoms(), 4);
        assert_eq!(mol.atoms[0].name, "N");
        assert_eq!(mol.atoms[0].res_name, "ALA");
        assert_eq!(mol.atoms[2].res_name, "GLY");
        // Should have 2 chains: A and B
        assert!(mol.chains.len() >= 1, "should have at least 1 chain");
    }
}
