//! Auto-detection of molecular file format from content.
//!
//! Reads the first few hundred bytes of a file and determines
//! whether it's PDB, SDF/MOL, XYZ, or mmCIF.
//!
//! Detection is based on content patterns, not file extensions.

use crate::types::FileFormat;

/// Detect the format of a molecular file from its content.
///
/// Reads up to the first 1KB to identify patterns:
/// - PDB: contains "ATOM  " or "HETATM" lines
/// - SDF/MOL: contains "V2000" or "V3000" count line, or "$$$$" delimiter
/// - XYZ: first non-empty line is an integer (atom count)
/// - mmCIF: starts with "data_" or contains "_atom_site."
pub fn detect_format(content: &str) -> Option<FileFormat> {
    let first_4k: String = content.chars().take(4096).collect();
    let lines: Vec<&str> = first_4k.lines().collect();

    // mmCIF FIRST: check for data_ header (takes priority over embedded ATOM-like text)
    let trimmed_start = first_4k.trim_start();
    if trimmed_start.starts_with("data_") || first_4k.contains("_atom_site.") {
        return Some(FileFormat::MmCif);
    }

    // PDB: look for ATOM/HETATM records in first 100 lines
    let atom_lines = lines.iter()
        .filter(|l| l.len() >= 6)
        .filter(|l| {
            let record = &l[..6];
            record.starts_with("ATOM") || record.starts_with("HETATM")
        })
        .count();

    if atom_lines >= 3 {
        return Some(FileFormat::Pdb);
    }

    // SDF/MOL: look for V2000/V3000 counts line or $$$$
    let has_delimiter = first_4k.contains("$$$$");
    let has_v2000 = first_4k.contains("V2000");
    let has_v3000 = first_4k.contains("V3000");
    if has_delimiter || has_v2000 || has_v3000 {
        return Some(FileFormat::Sdf);
    }

    // XYZ: first non-blank line is an integer (atom count)
    for line in &lines {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Ok(n) = trimmed.parse::<usize>() {
            if n > 0 && n < 10_000_000 {
                return Some(FileFormat::Xyz);
            }
        }
        break; // Only check first non-blank line
    }

    None
}

/// Detect format and also return the format name.
pub fn detect_format_name(content: &str) -> &'static str {
    match detect_format(content) {
        Some(FileFormat::Pdb) => "PDB",
        Some(FileFormat::Sdf) => "SDF/MOL",
        Some(FileFormat::Xyz) => "XYZ",
        Some(FileFormat::MmCif) => "mmCIF",
        Some(FileFormat::Mol) => "MOL",
        Some(FileFormat::Mol2) => "MOL2",
        None => "unknown",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detect_pdb() {
        let pdb = "HEADER\nATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N\nATOM      2  CA  ALA A   1      0.000   0.000   0.000  1.00  0.00           C\nATOM      3  C   ALA A   1      0.000   0.000   0.000  1.00  0.00           C\nEND\n";
        assert_eq!(detect_format(pdb), Some(FileFormat::Pdb));
    }

    #[test]
    fn detect_sdf() {
        let sdf = "Test\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n";
        assert_eq!(detect_format(sdf), Some(FileFormat::Sdf));
    }

    #[test]
    fn detect_xyz() {
        let xyz = "3\nWater molecule\nO  0.000  0.000  0.000\nH  0.957  0.000  0.000\nH -0.240  0.927  0.000\n";
        assert_eq!(detect_format(xyz), Some(FileFormat::Xyz));
    }

    #[test]
    fn detect_mmcif() {
        let cif = "data_1CRN\n#\n_atom_site.group_PDB\n";
        assert_eq!(detect_format(cif), Some(FileFormat::MmCif));
    }

    #[test]
    fn detect_unknown() {
        assert_eq!(detect_format("random text\nnot a molecule\n"), None);
    }

    #[test]
    fn detect_empty() {
        assert_eq!(detect_format(""), None);
    }
}
