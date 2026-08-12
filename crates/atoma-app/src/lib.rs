//! atoma desktop app — Tauri backend with JSON IPC.

use atoma_core::{parser::pdb::parse_pdb_str, analysis, detect_format};
use serde::Serialize;

#[derive(Serialize)]
struct MoleculeInfo {
    name: String, atoms: usize, residues: usize, chains: usize,
    mw: f64, format: String, sequence: String,
    bonds: usize, clashes: usize, rama_outliers: usize,
    ss_helix: usize, ss_strand: usize, ss_loop: usize, entropy: f64,
    ok: bool, error: String,
}

#[tauri::command]
fn open_file(content: String, filename: String) -> MoleculeInfo {
    let fmt = detect_format(&content).or_else(|| {
        // Fall back to extension if content detection fails
        let ext = std::path::Path::new(&filename).extension()?.to_str()?;
        match ext.to_lowercase().as_str() {
            "pdb" => Some(atoma_core::FileFormat::Pdb),
            "sdf" | "mol" => Some(atoma_core::FileFormat::Sdf),
            "xyz" => Some(atoma_core::FileFormat::Xyz),
            "cif" | "mmcif" => Some(atoma_core::FileFormat::MmCif),
            _ => None,
        }
    });
    let result = match fmt {
        Some(atoma_core::FileFormat::Pdb) => parse_pdb_str(&content, &filename)
            .map_err(|e| format!("{}", e)),
        Some(atoma_core::FileFormat::Sdf) |
        Some(atoma_core::FileFormat::Mol) => atoma_core::parse_sdf_str(&content)
            .map(|v| v.into_iter().next().unwrap())
            .map_err(|e| format!("{}", e)),
        Some(atoma_core::FileFormat::Xyz) => atoma_core::parse_xyz_str(&content)
            .map_err(|e| format!("{}", e)),
        Some(atoma_core::FileFormat::MmCif) => atoma_core::parse_mmcif_str(&content)
            .map_err(|e| format!("{}", e)),
        Some(atoma_core::FileFormat::Mol2) |
        None => Err("Unknown or unsupported format".into()),
    };

    match result {
        Ok(mol) => {
            let seqs = analysis::extract_sequences(&mol);
            let fasta = seqs.first().map(|(_, s)| s.clone()).unwrap_or_default();
            let name = mol.name.as_ref().cloned().unwrap_or(filename);
            MoleculeInfo {
                atoms: mol.n_atoms(), residues: mol.n_residues(), chains: mol.chains.len(),
                mw: mol.molecular_weight(), format: format!("{:?}", mol.source_format),
                sequence: fasta, name,
                bonds: analysis::detect_bonds(&mol).len(),
                clashes: analysis::detect_clashes(&mol).len(),
                rama_outliers: analysis::ramachandran(&mol).iter().filter(|r| !r.4).count(),
                ss_helix: analysis::secondary_structure(&mol).iter().filter(|s| s.2=='H').count(),
                ss_strand: analysis::secondary_structure(&mol).iter().filter(|s| s.2=='E').count(),
                ss_loop: analysis::secondary_structure(&mol).iter().filter(|s| s.2=='L').count(),
                entropy: analysis::structure_entropy(&mol),
                ok: true, error: String::new(),
            }
        }
        Err(e) => MoleculeInfo {
            name: filename, atoms: 0, residues: 0, chains: 0, mw: 0.0,
            format: "?".into(), sequence: String::new(),
            bonds: 0, clashes: 0, rama_outliers: 0,
            ss_helix: 0, ss_strand: 0, ss_loop: 0, entropy: 0.0,
            ok: false, error: e,
        }
    }
}

pub fn run() {
    tauri::Builder::default()
        .invoke_handler(tauri::generate_handler![open_file])
        .run(tauri::generate_context!())
        .expect("error");
}
