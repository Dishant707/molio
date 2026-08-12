use serde::{Deserialize, Serialize};

/// An individual atom in a molecular structure.
///
/// Coordinates are in Ångströms (Å), following the PDB convention.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Atom {
    /// Atom serial number (1-based, as in PDB)
    pub serial: u32,
    /// Atom name (e.g., "CA", "CB", "N")
    pub name: String,
    /// Alternate location indicator (' ' or 'A', 'B', ...)
    pub alt_loc: char,
    /// Residue name (e.g., "ALA", "GLY")
    pub res_name: String,
    /// Chain identifier
    pub chain_id: char,
    /// Residue sequence number
    pub res_seq: i32,
    /// Insertion code
    pub i_code: char,
    /// X coordinate in Å
    pub x: f64,
    /// Y coordinate in Å
    pub y: f64,
    /// Z coordinate in Å
    pub z: f64,
    /// Occupancy
    pub occupancy: f64,
    /// Temperature factor (B-factor)
    pub temp_factor: f64,
    /// Element symbol (e.g., "C", "N", "O", "S")
    pub element: String,
    /// Charge (if present)
    pub charge: Option<String>,
}

/// A residue (amino acid, nucleotide, ligand, etc.)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Residue {
    pub name: String,
    pub seq_num: i32,
    pub insertion_code: char,
    pub atoms: Vec<Atom>,
}

/// A chain in a macromolecular structure.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Chain {
    pub id: char,
    pub residues: Vec<Residue>,
}

/// A connection record (bond) between two atoms from PDB CONECT records.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Connection {
    pub atom1: (String, i32, char),  // (name, res_seq, chain)
    pub atom2: (String, i32, char),
}

/// Bond order / type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
    Unknown(u8),
}

/// A chemical bond between two atoms (1-based indices).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Bond {
    /// Index of first atom (1-based)
    pub atom1: u32,
    /// Index of second atom (1-based)
    pub atom2: u32,
    /// Bond order
    pub order: BondOrder,
}

/// Represents a parsed molecular structure.
///
/// Can be a macromolecule (protein/DNA) or a small molecule,
/// depending on the input format.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Molecule {
    /// Molecule name / title (from SDF header or PDB COMPND)
    pub name: Option<String>,
    /// PDB header / title lines
    pub header: Vec<String>,
    /// All atoms in the structure
    pub atoms: Vec<Atom>,
    /// Chains (populated for PDB/mmCIF macromolecules)
    pub chains: Vec<Chain>,
    /// Bond connectivity records (MOL/SDF bond block)
    pub bonds: Vec<Bond>,
    /// PDB CONECT connection records
    pub connections: Vec<Connection>,
    /// SDF properties (key-value pairs from > <tag> blocks)
    pub properties: std::collections::HashMap<String, String>,
    /// Original source format
    pub source_format: FileFormat,
}

/// Supported molecular file formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum FileFormat {
    Pdb,
    MmCif,
    Xyz,
    Sdf,
    Mol,
    Mol2,
}

impl Molecule {
    /// Create an empty molecule.
    pub fn new(format: FileFormat) -> Self {
        Molecule {
            name: None,
            header: Vec::new(),
            atoms: Vec::new(),
            chains: Vec::new(),
            bonds: Vec::new(),
            connections: Vec::new(),
            properties: std::collections::HashMap::new(),
            source_format: format,
        }
    }

    /// Total number of atoms.
    pub fn n_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// Total number of residues across all chains.
    pub fn n_residues(&self) -> usize {
        self.chains.iter().map(|c| c.residues.len()).sum()
    }

    /// Get the bounding box of all atom coordinates.
    pub fn bounding_box(&self) -> Option<([f64; 3], [f64; 3])> {
        if self.atoms.is_empty() {
            return None;
        }
        let mut min = [f64::MAX; 3];
        let mut max = [f64::MIN; 3];
        for atom in &self.atoms {
            min[0] = min[0].min(atom.x);
            min[1] = min[1].min(atom.y);
            min[2] = min[2].min(atom.z);
            max[0] = max[0].max(atom.x);
            max[1] = max[1].max(atom.y);
            max[2] = max[2].max(atom.z);
        }
        Some((min, max))
    }

    /// Get approximate molecular weight (sum of standard atomic weights).
    pub fn molecular_weight(&self) -> f64 {
        self.atoms.iter().map(|a| atomic_weight(&a.element)).sum()
    }
}

/// Standard atomic weights (g/mol) for common elements.
fn atomic_weight(element: &str) -> f64 {
    match element.trim() {
        "H" | "D" => 1.008,
        "C" => 12.011,
        "N" => 14.007,
        "O" => 16.000,
        "F" => 19.000,
        "P" => 30.974,
        "S" => 32.065,
        "CL" => 35.453,
        "MG" => 24.305,
        "CA" => 40.078,
        "ZN" => 65.380,
        "FE" => 55.845,
        "MN" => 54.938,
        "NA" => 22.990,
        "K" => 39.098,
        "SE" => 78.960,
        "BR" => 79.904,
        "I" => 126.904,
        _ => 12.0, // default to carbon weight for unknown
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_molecule() {
        let mol = Molecule::new(FileFormat::Pdb);
        assert_eq!(mol.n_atoms(), 0);
        assert_eq!(mol.n_residues(), 0);
        assert!(mol.bounding_box().is_none());
    }

    #[test]
    fn test_bounding_box() {
        let mut mol = Molecule::new(FileFormat::Pdb);
        mol.atoms.push(Atom {
            serial: 1, name: "CA".into(), alt_loc: ' ',
            res_name: "ALA".into(), chain_id: 'A', res_seq: 1,
            i_code: ' ', x: 0.0, y: 0.0, z: 0.0,
            occupancy: 1.0, temp_factor: 0.0,
            element: "C".into(), charge: None,
        });
        mol.atoms.push(Atom {
            serial: 2, name: "CA".into(), alt_loc: ' ',
            res_name: "ALA".into(), chain_id: 'A', res_seq: 2,
            i_code: ' ', x: 10.0, y: 20.0, z: 30.0,
            occupancy: 1.0, temp_factor: 0.0,
            element: "C".into(), charge: None,
        });

        let (min, max) = mol.bounding_box().unwrap();
        assert_eq!(min, [0.0, 0.0, 0.0]);
        assert_eq!(max, [10.0, 20.0, 30.0]);
    }
}
