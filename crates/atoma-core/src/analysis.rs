//! Molecular analysis utilities.

use crate::types::{Atom, Molecule};
use std::collections::HashMap;

/// Convert 3-letter amino acid code to 1-letter code.
fn aa3to1(code: &str) -> char {
    match code {
        "ALA" => 'A', "ARG" => 'R', "ASN" => 'N', "ASP" => 'D',
        "CYS" => 'C', "GLN" => 'Q', "GLU" => 'E', "GLY" => 'G',
        "HIS" => 'H', "ILE" => 'I', "LEU" => 'L', "LYS" => 'K',
        "MET" => 'M', "PHE" => 'F', "PRO" => 'P', "SER" => 'S',
        "THR" => 'T', "TRP" => 'W', "TYR" => 'Y', "VAL" => 'V',
        "ASX" => 'B', "GLX" => 'Z', "UNK" => 'X',
        // DNA/RNA
        "DA"  => 'A', "DC" => 'C', "DG" => 'G', "DT" => 'T',
        "A"   => 'A', "C"  => 'C', "G"  => 'G', "U"  => 'U',
        _ => 'X',
    }
}

/// Extract amino acid sequence from a molecule (FASTA format).
/// Returns one sequence per chain.
pub fn extract_sequences(mol: &Molecule) -> Vec<(char, String)> {
    let mut result = Vec::new();

    for chain in &mol.chains {
        let mut seq = String::new();
        let mut last_seq = -999;
        for residue in &chain.residues {
            // Avoid duplicates (same residue from alternate conformations)
            if residue.seq_num == last_seq && residue.insertion_code == ' ' {
                continue;
            }
            last_seq = residue.seq_num;
            seq.push(aa3to1(&residue.name));
        }
        if !seq.is_empty() {
            result.push((chain.id, seq));
        }
    }

    // If no chains, try to extract from flat atom list
    if result.is_empty() && !mol.atoms.is_empty() {
        let mut seq = String::new();
        let mut last_res = String::new();
        let mut last_seq = -999;
        for atom in &mol.atoms {
            if atom.res_seq == last_seq && atom.res_name == last_res {
                continue;
            }
            if !atom.res_name.is_empty() {
                seq.push(aa3to1(&atom.res_name));
                last_res = atom.res_name.clone();
                last_seq = atom.res_seq;
            }
        }
        if !seq.is_empty() {
            result.push(('A', seq));
        }
    }

    result
}

/// Format sequences as FASTA string.
pub fn to_fasta(sequences: &[(char, String)], header: &str) -> String {
    let mut fasta = String::new();
    for (chain_id, seq) in sequences {
        fasta.push_str(&format!(">{}:{}|atoma\n", header, chain_id));
        // Wrap at 60 chars per line (FASTA convention)
        for chunk in seq.as_bytes().chunks(60) {
            fasta.push_str(std::str::from_utf8(chunk).unwrap());
            fasta.push('\n');
        }
    }
    fasta
}

// ─── Geometry-Based Bond Detection ───────────────────────────────

/// Covalent bond distance thresholds (Å) for element pairs.
fn bond_range(e1: &str, e2: &str) -> Option<(f64, f64)> {
    let key = if e1 < e2 { (e1, e2) } else { (e2, e1) };
    match key {
        ("C", "C") => Some((1.2, 1.7)), ("C", "N") => Some((1.2, 1.55)),
        ("C", "O") => Some((1.15, 1.55)), ("C", "S") => Some((1.6, 1.9)),
        ("H", "C") => Some((0.9, 1.2)), ("H", "N") => Some((0.8, 1.15)),
        ("H", "O") => Some((0.8, 1.1)), ("N", "N") => Some((1.1, 1.5)),
        ("N", "O") => Some((1.1, 1.5)), ("O", "O") => Some((1.1, 1.55)),
        ("O", "S") => Some((1.4, 1.8)), ("S", "S") => Some((1.85, 2.2)),
        ("C", "P") => Some((1.6, 1.95)), ("O", "P") => Some((1.4, 1.75)),
        _ => None,
    }
}

/// Detect bonds from atom distances.
pub fn detect_bonds(mol: &Molecule) -> Vec<(usize, usize, f64)> {
    let mut bonds = Vec::new();
    let n = mol.n_atoms();
    for i in 0..n {
        for j in (i + 1)..n {
            let a = &mol.atoms[i]; let b = &mol.atoms[j];
            let dx = a.x - b.x; let dy = a.y - b.y; let dz = a.z - b.z;
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            if let Some((min_d, max_d)) = bond_range(&a.element, &b.element) {
                if dist >= min_d && dist <= max_d { bonds.push((i, j, dist)); }
            }
        }
    }
    bonds
}

// ─── Steric Clash Detection ──────────────────────────────────────

fn vdw_radius(element: &str) -> f64 {
    match element {
        "H" => 1.20, "C" => 1.70, "N" => 1.55, "O" => 1.52,
        "F" => 1.47, "P" => 1.80, "S" => 1.80, "Cl" => 1.75, "Br" => 1.85,
        "Fe" => 2.00, "Zn" => 1.90, "Mg" => 1.73, "Ca" => 2.00,
        _ => 1.70,
    }
}

/// Detect steric clashes (atoms closer than 0.8 × sum of vdW radii).
pub fn detect_clashes(mol: &Molecule) -> Vec<(usize, usize, f64, f64)> {
    let mut clashes = Vec::new();
    let n = mol.n_atoms();
    let radii: Vec<f64> = mol.atoms.iter().map(|a| vdw_radius(&a.element)).collect();
    for i in 0..n {
        for j in (i + 1)..n {
            if (j as i32 - i as i32).abs() <= 3 { continue; }
            let a = &mol.atoms[i]; let b = &mol.atoms[j];
            let dx = a.x - b.x; let dy = a.y - b.y; let dz = a.z - b.z;
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            let threshold = 0.8 * (radii[i] + radii[j]);
            if dist < threshold && dist > 0.1 { clashes.push((i, j, dist, threshold)); }
        }
    }
    clashes
}

// ─── Ramachandran Analysis ───────────────────────────────────────

fn dihedral(p1: &[f64; 3], p2: &[f64; 3], p3: &[f64; 3], p4: &[f64; 3]) -> f64 {
    let b1 = [p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]];
    let b2 = [p3[0]-p2[0], p3[1]-p2[1], p3[2]-p2[2]];
    let b3 = [p4[0]-p3[0], p4[1]-p3[1], p4[2]-p3[2]];
    let n1 = [b1[1]*b2[2]-b1[2]*b2[1], b1[2]*b2[0]-b1[0]*b2[2], b1[0]*b2[1]-b1[1]*b2[0]];
    let n2 = [b2[1]*b3[2]-b2[2]*b3[1], b2[2]*b3[0]-b2[0]*b3[2], b2[0]*b3[1]-b2[1]*b3[0]];
    let n1n = (n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2]).sqrt();
    let n2n = (n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2]).sqrt();
    if n1n < 1e-10 || n2n < 1e-10 { return 0.0; }
    let cos = ((n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2])/(n1n*n2n)).clamp(-1.0, 1.0);
    let sign = if b2[0]*(n1[1]*n2[2]-n1[2]*n2[1])+b2[1]*(n1[2]*n2[0]-n1[0]*n2[2])+b2[2]*(n1[0]*n2[1]-n1[1]*n2[0]) >= 0.0 {1.0} else {-1.0};
    sign * cos.acos().to_degrees()
}

pub fn ramachandran(mol: &Molecule) -> Vec<(String, i32, f64, f64, bool)> {
    let mut results = Vec::new();
    for chain in &mol.chains {
        for window in chain.residues.windows(3) {
            let rp=&window[0]; let rc=&window[1]; let rn=&window[2];
            let cp=rp.atoms.iter().find(|a| a.name=="C");
            let n=rc.atoms.iter().find(|a| a.name=="N");
            let ca=rc.atoms.iter().find(|a| a.name=="CA");
            let c=rc.atoms.iter().find(|a| a.name=="C");
            let nn=rn.atoms.iter().find(|a| a.name=="N");
            if let (Some(cp),Some(n),Some(ca),Some(c),Some(nn))=(cp,n,ca,c,nn) {
                let phi=dihedral(&[cp.x,cp.y,cp.z],&[n.x,n.y,n.z],&[ca.x,ca.y,ca.z],&[c.x,c.y,c.z]);
                let psi=dihedral(&[n.x,n.y,n.z],&[ca.x,ca.y,ca.z],&[c.x,c.y,c.z],&[nn.x,nn.y,nn.z]);
                let allowed = (phi > -100. && phi < -30. && psi > -80. && psi < -10.)
                    || (phi > -180. && phi < -50. && psi > 80.)
                    || (phi > 20. && phi < 100. && psi > -20. && psi < 80.);
                results.push((rc.name.clone(),rc.seq_num,phi,psi,allowed));
            }
        }
    }
    results
}

// ─── Secondary Structure ─────────────────────────────────────────

pub fn secondary_structure(mol: &Molecule) -> Vec<(String, i32, char)> {
    let mut result = Vec::new();
    for chain in &mol.chains {
        let cas: Vec<&crate::types::Atom> = chain.residues.iter()
            .filter_map(|r| r.atoms.iter().find(|a| a.name=="CA")).collect();
        if cas.len() < 5 { continue; }
        for i in 0..cas.len() {
            let ss = if i+4 < cas.len() {
                let ca=&cas[i]; let ca4=&cas[i+4];
                let d=(ca.x-ca4.x).powi(2)+(ca.y-ca4.y).powi(2)+(ca.z-ca4.z).powi(2);
                if d.sqrt() < 6.5 {'H'} else {'L'}
            } else {'L'};
            result.push((chain.residues[i].name.clone(), chain.residues[i].seq_num, ss));
        }
    }
    result
}

// ─── Combined Analysis ───────────────────────────────────────────

pub fn analyze(mol: &Molecule) -> String {
    let mut r = String::new();
    r.push_str("╔══════════════════════════════════╗\n");
    r.push_str("║  ⚛️  atoma — Structural Analysis   ║\n");
    r.push_str("╠══════════════════════════════════╣\n");
    r.push_str(&format!("║  Atoms:     {:>22} ║\n", mol.n_atoms()));
    if !mol.chains.is_empty() {
        r.push_str(&format!("║  Residues:  {:>22} ║\n", mol.n_residues()));
        r.push_str(&format!("║  Chains:    {:>22} ║\n", mol.chains.len()));
    }
    r.push_str("╠══════════════════════════════════╣\n");
    let h = structure_entropy(mol);
    r.push_str(&format!("║  Structure entropy: {:>13.2} bits ║\n", h));
    let bonds = detect_bonds(mol);
    r.push_str(&format!("║  Geometry bonds:  {:>13} ║\n", bonds.len()));
    let clashes = detect_clashes(mol);
    if clashes.is_empty() { r.push_str("║  ✅ No steric clashes             ║\n"); }
    else { r.push_str(&format!("║  ⚠️  {} steric clashes             ║\n", clashes.len())); }
    let rama = ramachandran(mol);
    if !rama.is_empty() {
        let outliers = rama.iter().filter(|r| !r.4).count();
        if outliers == 0 { r.push_str("║  ✅ Ramachandran: 0 outliers      ║\n"); }
        else { r.push_str(&format!("║  ⚠️  Rama outliers: {:>13} ║\n", outliers)); }
    }
    let ss = secondary_structure(mol);
    if !ss.is_empty() {
        let h = ss.iter().filter(|s| s.2=='H').count();
        let e = ss.iter().filter(|s| s.2=='E').count();
        r.push_str(&format!("║  SS: H={} E={} L={}          ║\n", h, e, ss.len()-h-e));
    }
    r.push_str("╚══════════════════════════════════╝\n");
    r
}

// ─── Information-Theoretic Analysis (Shannon, 1948) ─────────────

/// Shannon entropy of a distribution: H = -∑ p_i · log₂(p_i)
/// Measures disorder/uncertainty in bits.
///
/// Reference: Shannon, C.E. (1948). "A Mathematical Theory of Communication."
/// Bell System Technical Journal, 27(3), 379-423.
///
/// Application to proteins: Chanda et al. (2020). "Information Theory in
/// Computational Biology." Entropy, 22(6), 627.
fn shannon_h(counts: &[f64]) -> f64 {
    let total: f64 = counts.iter().sum();
    if total == 0.0 { return 0.0; }
    counts.iter()
        .filter(|&&c| c > 0.0)
        .map(|&c| {
            let p = c / total;
            -p * p.log2()
        })
        .sum()
}

/// Per-residue Shannon entropy from atom type distribution.
/// High entropy (>2.0) = flexible/disordered residue.
/// Low entropy (<1.0) = rigid/well-ordered residue.
pub fn residue_entropy(mol: &Molecule) -> Vec<(String, i32, f64)> {
    let mut results = Vec::new();
    for chain in &mol.chains {
        for residue in &chain.residues {
            let mut atom_counts: HashMap<&str, f64> = HashMap::new();
            for atom in &residue.atoms {
                *atom_counts.entry(&atom.name).or_insert(0.0) += 1.0;
            }
            let counts: Vec<f64> = atom_counts.values().copied().collect();
            let h = shannon_h(&counts);
            results.push((residue.name.clone(), residue.seq_num, h));
        }
    }
    results
}

/// Information density of the entire structure.
/// High values suggest diverse atom types (complex molecule).
/// Low values suggest repetitive structure.
pub fn structure_entropy(mol: &Molecule) -> f64 {
    let mut element_counts: HashMap<&str, f64> = HashMap::new();
    for atom in &mol.atoms {
        *element_counts.entry(&atom.element).or_insert(0.0) += 1.0;
    }
    let counts: Vec<f64> = element_counts.values().copied().collect();
    shannon_h(&counts)
}

/// B-factor entropy — measures thermal disorder from crystallographic data.
/// Reference: Strait & Dewey (1996). "Shannon Information Entropy of
/// Protein Sequences." Biophysical Journal, 71(1), 148-155.
pub fn bfactor_entropy(mol: &Molecule) -> f64 {
    let bfactors: Vec<f64> = mol.atoms.iter()
        .filter(|a| a.temp_factor > 0.0)
        .map(|a| a.temp_factor)
        .collect();
    if bfactors.is_empty() { return 0.0; }
    // Bin B-factors into 10 bins for entropy calculation
    let max_b = bfactors.iter().cloned().fold(0.0_f64, f64::max);
    let min_b = bfactors.iter().cloned().fold(f64::MAX, f64::min);
    let range = (max_b - min_b).max(1.0);
    let mut bins = [0.0_f64; 10];
    for b in &bfactors {
        let idx = ((b - min_b) / range * 9.999) as usize;
        bins[idx.min(9)] += 1.0;
    }
    shannon_h(&bins)
}

// ─── Format Writers ──────────────────────────────────────────────

/// Write molecule as PDB format.
pub fn write_pdb(mol: &Molecule) -> String {
    let mut pdb = String::new();
    for atom in &mol.atoms {
        pdb.push_str(&format!(
            "ATOM  {:>5} {:^4}{}{:3} {}{:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}  \n",
            atom.serial,
            atom.name,
            if atom.alt_loc == ' ' { ' ' } else { atom.alt_loc },
            atom.res_name,
            atom.chain_id,
            atom.res_seq,
            atom.x, atom.y, atom.z,
            atom.occupancy, atom.temp_factor,
            atom.element,
        ));
    }
    pdb.push_str("END\n");
    pdb
}

/// Write molecule as XYZ format.
pub fn write_xyz(mol: &Molecule, comment: &str) -> String {
    let mut xyz = String::new();
    xyz.push_str(&format!("{}\n", mol.n_atoms()));
    xyz.push_str(&format!("{}\n", comment));
    for atom in &mol.atoms {
        xyz.push_str(&format!(
            "{:<3} {:>12.6} {:>12.6} {:>12.6}\n",
            atom.element, atom.x, atom.y, atom.z
        ));
    }
    xyz
}

/// Write molecule as SDF/MOL format.
pub fn write_sdf(mol: &Molecule) -> String {
    let n_atoms = mol.n_atoms();
    let n_bonds = mol.bonds.len();

    let mut sdf = String::new();
    let name = mol.name.as_deref().unwrap_or("molio");
    sdf.push_str(&format!("{}\n", name));
    sdf.push_str("  molio\n\n");
    sdf.push_str(&format!("{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n", n_atoms, n_bonds));

    for atom in &mol.atoms {
        sdf.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
            atom.x, atom.y, atom.z, atom.element
        ));
    }

    for bond in &mol.bonds {
        let bond_type = match bond.order {
            crate::types::BondOrder::Single => 1,
            crate::types::BondOrder::Double => 2,
            crate::types::BondOrder::Triple => 3,
            crate::types::BondOrder::Aromatic => 4,
            crate::types::BondOrder::Unknown(n) => n as i32,
        };
        sdf.push_str(&format!("{:>3}{:>3}{:>3}  0  0  0  0\n", bond.atom1, bond.atom2, bond_type));
    }

    sdf.push_str("M  END\n");
    for (key, val) in &mol.properties {
        sdf.push_str(&format!("> <{}>\n{}\n\n", key, val));
    }
    sdf.push_str("$$$$\n");
    sdf
}

/// Convert between formats.
pub fn convert_format(mol: &Molecule, target: &str) -> Option<String> {
    match target.to_lowercase().as_str() {
        "pdb" => Some(write_pdb(mol)),
        "xyz" => Some(write_xyz(mol, "Converted by atoma")),
        "sdf" | "mol" => Some(write_sdf(mol)),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::pdb::parse_pdb_str;

    #[test]
    fn test_extract_sequence_1crn() {
        let pdb = include_str!("../../../test_data/pdb/1crn.pdb");
        let mol = parse_pdb_str(pdb, "t.pdb").unwrap();
        let seqs = extract_sequences(&mol);
        assert_eq!(seqs.len(), 1);
        assert_eq!(seqs[0].0, 'A');
        // 1CRN: TTCCPSIVAR... (10 residues)
        assert!(seqs[0].1.starts_with("TT"));
        assert_eq!(seqs[0].1.len(), 10);
    }

    #[test]
    fn test_to_fasta() {
        let seqs = vec![('A', "TTCCPSIVAR".to_string())];
        let fasta = to_fasta(&seqs, "1crn");
        assert!(fasta.starts_with(">1crn:A|atoma\n"));
        assert!(fasta.contains("TTCCPSIVAR"));
    }

    #[test]
    fn test_write_pdb_roundtrip() {
        let pdb = include_str!("../../../test_data/pdb/1crn.pdb");
        let mol = parse_pdb_str(pdb, "t.pdb").unwrap();
        let written = write_pdb(&mol);
        let mol2 = parse_pdb_str(&written, "r.pdb").unwrap();
        assert_eq!(mol.n_atoms(), mol2.n_atoms());
        // First atom should match
        assert!((mol.atoms[0].x - mol2.atoms[0].x).abs() < 0.01);
    }

    #[test]
    fn test_write_xyz() {
        let pdb = include_str!("../../../test_data/pdb/1crn.pdb");
        let mol = parse_pdb_str(pdb, "t.pdb").unwrap();
        let xyz = write_xyz(&mol, "test");
        let first_line: String = xyz.lines().next().unwrap().to_string();
        assert_eq!(first_line, mol.n_atoms().to_string());
    }
}
