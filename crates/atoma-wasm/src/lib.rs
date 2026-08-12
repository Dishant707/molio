//! atoma WebAssembly bindings.
//! Compile with: wasm-pack build --target web

use atoma_core::{detect_format, parser};
use serde::Serialize;
use wasm_bindgen::prelude::*;

#[derive(Serialize)]
pub struct MoleculeInfo {
    pub name: String,
    pub atoms: usize,
    pub residues: usize,
    pub chains: usize,
    pub mw: f64,
    pub format: String,
    pub sequence: String,
    pub bonds: usize,
    pub clashes: usize,
    pub rama_outliers: usize,
    pub ss_helix: usize,
    pub ss_strand: usize,
    pub ss_loop: usize,
    pub entropy: f64,
    pub ok: bool,
    pub error: String,
}

#[wasm_bindgen]
pub fn parse_file(content: &str, filename: &str) -> JsValue {
    let fmt = detect_format(content).or_else(|| {
        let ext = std::path::Path::new(filename)
            .extension()
            .and_then(|e| e.to_str())?;
        match ext.to_lowercase().as_str() {
            "pdb" => Some(atoma_core::FileFormat::Pdb),
            "sdf" | "mol" => Some(atoma_core::FileFormat::Sdf),
            "xyz" => Some(atoma_core::FileFormat::Xyz),
            "cif" | "mmcif" => Some(atoma_core::FileFormat::MmCif),
            _ => None,
        }
    });

    let result = match fmt {
        Some(atoma_core::FileFormat::Pdb) => {
            parser::pdb::parse_pdb_str(content, filename).map_err(|e| format!("{e}"))
        }
        Some(atoma_core::FileFormat::Sdf) | Some(atoma_core::FileFormat::Mol) => {
            atoma_core::parse_sdf_str(content)
                .map(|v| v.into_iter().next().unwrap())
                .map_err(|e| format!("{e}"))
        }
        Some(atoma_core::FileFormat::Xyz) => {
            atoma_core::parse_xyz_str(content).map_err(|e| format!("{e}"))
        }
        Some(atoma_core::FileFormat::MmCif) => {
            atoma_core::parse_mmcif_str(content).map_err(|e| format!("{e}"))
        }
        _ => Err("Unknown format".into()),
    };

    let info = match result {
        Ok(mol) => {
            let seqs = atoma_core::analysis::extract_sequences(&mol);
            let fasta = seqs.first().map(|(_, s)| s.clone()).unwrap_or_default();
            let name = mol.name.as_ref().cloned().unwrap_or_else(|| filename.to_string());
            MoleculeInfo {
                atoms: mol.n_atoms(),
                residues: mol.n_residues(),
                chains: mol.chains.len(),
                mw: mol.molecular_weight(),
                format: format!("{:?}", mol.source_format),
                sequence: fasta,
                name,
                bonds: atoma_core::analysis::detect_bonds(&mol).len(),
                clashes: atoma_core::analysis::detect_clashes(&mol).len(),
                rama_outliers: atoma_core::analysis::ramachandran(&mol)
                    .iter()
                    .filter(|r| !r.4)
                    .count(),
                ss_helix: atoma_core::analysis::secondary_structure(&mol)
                    .iter()
                    .filter(|s| s.2 == 'H')
                    .count(),
                ss_strand: atoma_core::analysis::secondary_structure(&mol)
                    .iter()
                    .filter(|s| s.2 == 'E')
                    .count(),
                ss_loop: atoma_core::analysis::secondary_structure(&mol)
                    .iter()
                    .filter(|s| s.2 == 'L')
                    .count(),
                entropy: atoma_core::analysis::structure_entropy(&mol),
                ok: true,
                error: String::new(),
            }
        }
        Err(e) => MoleculeInfo {
            name: filename.to_string(),
            atoms: 0,
            residues: 0,
            chains: 0,
            mw: 0.0,
            format: "?".into(),
            sequence: String::new(),
            bonds: 0,
            clashes: 0,
            rama_outliers: 0,
            ss_helix: 0,
            ss_strand: 0,
            ss_loop: 0,
            entropy: 0.0,
            ok: false,
            error: e,
        },
    };

    serde_wasm_bindgen::to_value(&info).unwrap_or(JsValue::NULL)
}
