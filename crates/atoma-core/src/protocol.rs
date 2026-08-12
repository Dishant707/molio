//! Binary communication protocol for molio frontend-backend bridge.
//!
//! Replaces JSON serialization with a compact binary format.
//! Based on zero-copy principles from:
//! - Lemire & Langdale (2019). "Parsing Gigabytes of JSON per Second." VLDB.
//! - Chen et al. (2024). "Lite2: Zero-Copy Serialization." Computers, 13(4):89.

/// Binary header: magic + version + field count
const MAGIC: u32 = 0x4D4F4C49; // "ATOM"
const VERSION: u16 = 1;

/// Write a molecule analysis result as compact binary.
/// Format: [4B magic][2B version][2B field_count][fields...]
/// Each field: [1B type][2B len][len bytes data]
pub fn write_analysis_binary(
    name: &str, atoms: u32, residues: u32, chains: u32,
    mw: f64, bonds: u32, clashes: u32, rama_outliers: u32,
    ss_helix: u32, ss_strand: u32, ss_loop: u32,
    sequence: &str, format_name: &str,
    structure_entropy: f64,
) -> Vec<u8> {
    let mut buf = Vec::with_capacity(256);

    // Header
    buf.extend_from_slice(&MAGIC.to_le_bytes());
    buf.extend_from_slice(&VERSION.to_le_bytes());

    // Field count (14 fields)
    let field_count: u16 = 14;
    buf.extend_from_slice(&field_count.to_le_bytes());

    // Helper: write field as [type][len][data]
    macro_rules! write_field {
        ($buf:expr, $type:expr, $data:expr) => {{
            $buf.push($type);
            let len = $data.len() as u16;
            $buf.extend_from_slice(&len.to_le_bytes());
            $buf.extend_from_slice($data);
        }};
    }

    write_field!(buf, 0x01, &atoms.to_le_bytes());       // u32
    write_field!(buf, 0x01, &residues.to_le_bytes());
    write_field!(buf, 0x01, &chains.to_le_bytes());
    write_field!(buf, 0x02, &mw.to_le_bytes());          // f64
    write_field!(buf, 0x01, &bonds.to_le_bytes());
    write_field!(buf, 0x01, &clashes.to_le_bytes());
    write_field!(buf, 0x01, &rama_outliers.to_le_bytes());
    write_field!(buf, 0x01, &ss_helix.to_le_bytes());
    write_field!(buf, 0x01, &ss_strand.to_le_bytes());
    write_field!(buf, 0x01, &ss_loop.to_le_bytes());
    write_field!(buf, 0x03, name.as_bytes());            // string
    write_field!(buf, 0x03, sequence.as_bytes());
    write_field!(buf, 0x03, format_name.as_bytes());
    write_field!(buf, 0x02, &structure_entropy.to_le_bytes());

    buf
}

/// Read a u32 from little-endian bytes.
fn read_u32(data: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes([data[offset], data[offset+1], data[offset+2], data[offset+3]])
}

/// Read a f64 from little-endian bytes.
fn read_f64(data: &[u8], offset: usize) -> f64 {
    f64::from_le_bytes([
        data[offset], data[offset+1], data[offset+2], data[offset+3],
        data[offset+4], data[offset+5], data[offset+6], data[offset+7],
    ])
}

/// Read a string field.
fn read_string(data: &[u8], offset: &mut usize) -> String {
    let len = u16::from_le_bytes([data[*offset], data[*offset+1]]) as usize;
    *offset += 2;
    let s = String::from_utf8_lossy(&data[*offset..*offset+len]).to_string();
    *offset += len;
    s
}

/// Parse binary analysis result. Returns parsed fields.
pub fn read_analysis_binary(data: &[u8]) -> Option<(
    u32, u32, u32, f64, u32, u32, u32, u32, u32, u32, String, String, String, f64
)> {
    if data.len() < 8 { return None; }

    let magic = read_u32(data, 0);
    if magic != MAGIC { return None; }

    let _version = u16::from_le_bytes([data[4], data[5]]);
    let field_count = u16::from_le_bytes([data[6], data[7]]);
    if field_count != 14 { return None; }

    let mut offset = 8;

    macro_rules! read_u32_field {
        () => {{
            let _type = data[offset]; offset += 1;
            let _len = u16::from_le_bytes([data[offset], data[offset+1]]); offset += 2;
            let val = read_u32(data, offset); offset += 4;
            val
        }};
    }

    macro_rules! read_f64_field {
        () => {{
            let _type = data[offset]; offset += 1;
            let _len = u16::from_le_bytes([data[offset], data[offset+1]]); offset += 2;
            let val = read_f64(data, offset); offset += 8;
            val
        }};
    }

    macro_rules! read_str_field {
        () => {{
            let _type = data[offset]; offset += 1;
            read_string(data, &mut offset)
        }};
    }

    let atoms = read_u32_field!();
    let residues = read_u32_field!();
    let chains = read_u32_field!();
    let mw = read_f64_field!();
    let bonds = read_u32_field!();
    let clashes = read_u32_field!();
    let rama_outliers = read_u32_field!();
    let ss_helix = read_u32_field!();
    let ss_strand = read_u32_field!();
    let ss_loop = read_u32_field!();
    let name = read_str_field!();
    let sequence = read_str_field!();
    let format_name = read_str_field!();
    let entropy = read_f64_field!();

    Some((atoms, residues, chains, mw, bonds, clashes, rama_outliers,
          ss_helix, ss_strand, ss_loop, name, sequence, format_name, entropy))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_binary_roundtrip() {
        let data = write_analysis_binary(
            "1crn", 71, 10, 1, 974.7, 71, 0, 0, 0, 0, 10,
            "TTCCPSIVRS", "Pdb", 1.85,
        );

        let parsed = read_analysis_binary(&data).unwrap();
        assert_eq!(parsed.0, 71);   // atoms
        assert_eq!(parsed.1, 10);   // residues
        assert!((parsed.3 - 974.7).abs() < 0.1); // mw
        assert_eq!(parsed.10, "1crn"); // name
        assert_eq!(parsed.11, "TTCCPSIVRS"); // sequence
    }
}
