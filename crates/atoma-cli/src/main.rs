use clap::{Parser, Subcommand};
use std::path::Path;
use std::time::Instant;

use atoma_core::{parse_pdb, parse_sdf, parse_xyz, parse_mmcif, detect_format, FileFormat};
use atoma_core::analysis;

/// Detect file format from content, falling back to extension.
fn detect(path: &str) -> FileFormat {
    // Try content-based detection first
    if let Ok(content) = std::fs::read_to_string(path) {
        if let Some(fmt) = detect_format(&content) {
            return fmt;
        }
    }
    // Fall back to extension
    let path = Path::new(path);
    match path.extension().and_then(|e| e.to_str()) {
        Some("pdb") => FileFormat::Pdb,
        Some("sdf") | Some("mol") => FileFormat::Sdf,
        Some("xyz") => FileFormat::Xyz,
        Some("cif") => FileFormat::MmCif,
        _ => FileFormat::Pdb,
    }
}

#[derive(Parser)]
#[command(name = "atoma", version, about = "⚛️ High-performance molecular file I/O")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Parse and display molecule information
    View {
        /// Path to the molecular file
        file: String,
    },
    /// Run benchmark against the file
    Bench {
        /// Path to the molecular file
        file: String,
        /// Number of iterations for warm-up
        #[arg(short, long, default_value = "10")]
        warmup: u32,
        /// Number of iterations to measure
        #[arg(short, long, default_value = "100")]
        iterations: u32,
    },
    /// Print molecule statistics
    Stats {
        /// Path to the molecular file
        file: String,
    },
    /// Extract amino acid sequence as FASTA
    ExtractSequence {
        /// Path to the molecular file
        file: String,
    },
    /// Convert between file formats
    Convert {
        /// Path to the input file
        file: String,
        /// Target format (pdb, xyz, sdf)
        #[arg(short, long)]
        to: String,
        /// Output file path (stdout if not specified)
        #[arg(short, long)]
        output: Option<String>,
    },
    /// Run structural analysis (bonds, clashes, Ramachandran, SS)
    Analyze {
        /// Path to the molecular file
        file: String,
    },
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::View { file } => cmd_view(&file),
        Commands::Bench { file, warmup, iterations } => cmd_bench(&file, warmup, iterations),
        Commands::Stats { file } => cmd_stats(&file),
        Commands::ExtractSequence { file } => cmd_extract_sequence(&file),
        Commands::Convert { file, to, output } => cmd_convert(&file, &to, output.as_deref()),
        Commands::Analyze { file } => cmd_analyze(&file),
    }
}

fn cmd_view(path: &str) -> anyhow::Result<()> {
    let fmt = detect(path);
    let start = Instant::now();

    match fmt {
        FileFormat::MmCif => {
            let mol = parse_mmcif(path)?;
            let elapsed = start.elapsed();
            print_molecule(&mol, 0, 1);
            println!("\n  Parse time: {:?}", elapsed);
        }
        FileFormat::Xyz => {
            let mol = parse_xyz(path)?;
            let elapsed = start.elapsed();
            print_molecule(&mol, 0, 1);
            println!("\n  Parse time: {:?}", elapsed);
        }
        FileFormat::Sdf | FileFormat::Mol => {
            let mols = parse_sdf(path)?;
            let elapsed = start.elapsed();
            for (i, mol) in mols.iter().enumerate() {
                print_molecule(mol, i + 1, mols.len());
            }
            println!("\n  Total: {} molecules in {:?}", mols.len(), elapsed);
        }
        _ => {
            let mol = parse_pdb(path)?;
            let elapsed = start.elapsed();
            print_molecule(&mol, 0, 1);
            println!("\n  Parse time: {:?}", elapsed);
        }
    }
    Ok(())
}

fn print_molecule(mol: &atoma_core::Molecule, idx: usize, total: usize) {
    if total > 1 {
        println!("\n  ═══ Molecule {}/{} ═══", idx, total);
    }
    println!("╔══════════════════════════════════╗");
    println!("║  ⚛️  atoma — Molecular Viewer     ║");
    println!("╠══════════════════════════════════╣");
    if let Some(ref name) = mol.name {
        println!("║ Name:     {:22} ║", truncate(name, 22));
    }
    println!("║ Format:   {:22} ║", format!("{:?}", mol.source_format));
    println!("║ Atoms:    {:>22} ║", mol.n_atoms());
    println!("║ Bonds:    {:>22} ║", mol.bonds.len());
    if !mol.chains.is_empty() {
        println!("║ Residues: {:>22} ║", mol.n_residues());
        println!("║ Chains:   {:>22} ║", mol.chains.len());
    }
    println!("║ MW:       {:>19.1} Da ║", mol.molecular_weight());
    println!("╚══════════════════════════════════╝");

    if !mol.properties.is_empty() {
        println!("\n  Properties:");
        for (k, v) in &mol.properties {
            println!("    {}: {}", k, v);
        }
    }

    if let Some((min, max)) = mol.bounding_box() {
        println!("\n  Bounding box:");
        println!("    min: [{:.3}, {:.3}, {:.3}]", min[0], min[1], min[2]);
        println!("    max: [{:.3}, {:.3}, {:.3}]", max[0], max[1], max[2]);
    }

    if !mol.chains.is_empty() {
        println!("\n  Chains:");
        for chain in &mol.chains {
            println!("    Chain {}: {} residues, {} atoms",
                chain.id,
                chain.residues.len(),
                chain.residues.iter().map(|r| r.atoms.len()).sum::<usize>()
            );
        }
    }
}

fn cmd_bench(path: &str, warmup: u32, iterations: u32) -> anyhow::Result<()> {
    let fmt = detect(path);
    println!("⚡ atoma benchmark");
    println!("   File:       {path}");
    println!("   Format:     {:?}", fmt);
    println!("   Warmup:     {warmup} iterations");
    println!("   Iterations: {iterations} iterations");
    println!();

    // Read file once for fairness
    let content = std::fs::read_to_string(path)?;
    let size_mb = content.len() as f64 / (1024.0 * 1024.0);

    match fmt {
        FileFormat::MmCif => {
            for _ in 0..warmup {
                let _ = atoma_core::parser::mmcif::parse_mmcif_str(&content);
            }
            let start = Instant::now();
            for _ in 0..iterations {
                let _ = atoma_core::parser::mmcif::parse_mmcif_str(&content);
            }
            let total = start.elapsed();
            let avg = total / iterations;
            let mol = atoma_core::parser::mmcif::parse_mmcif_str(&content)?;

            println!("  ┌─────────────────────────────────┐");
            println!("  │  📊 Results                     │");
            println!("  ├─────────────────────────────────┤");
            println!("  │  Atoms:     {:>18}      │", mol.n_atoms());
            println!("  │  File size: {:>15.2} MB      │", size_mb);
            println!("  │  Total:     {:>18.2?}      │", total);
            println!("  │  Average:   {:>18.2?}      │", avg);
            println!("  │  Throughput:{:>15.1} atoms/ms  │",
                mol.n_atoms() as f64 / avg.as_secs_f64() / 1000.0
            );
            println!("  └─────────────────────────────────┘");
        }
        FileFormat::Xyz => {
            for _ in 0..warmup {
                let _ = atoma_core::parser::xyz::parse_xyz_str(&content);
            }
            let start = Instant::now();
            for _ in 0..iterations {
                let _ = atoma_core::parser::xyz::parse_xyz_str(&content);
            }
            let total = start.elapsed();
            let avg = total / iterations;
            let mol = atoma_core::parser::xyz::parse_xyz_str(&content)?;

            println!("  ┌─────────────────────────────────┐");
            println!("  │  📊 Results                     │");
            println!("  ├─────────────────────────────────┤");
            println!("  │  Atoms:     {:>18}      │", mol.n_atoms());
            println!("  │  File size: {:>15.2} MB      │", size_mb);
            println!("  │  Total:     {:>18.2?}      │", total);
            println!("  │  Average:   {:>18.2?}      │", avg);
            println!("  │  Throughput:{:>15.1} atoms/ms  │",
                mol.n_atoms() as f64 / avg.as_secs_f64() / 1000.0
            );
            println!("  └─────────────────────────────────┘");
        }
        FileFormat::Sdf | FileFormat::Mol => {
            // Warmup
            for _ in 0..warmup {
                let _ = atoma_core::parser::sdf::parse_sdf_str(&content);
            }
            // Benchmark
            let start = Instant::now();
            for _ in 0..iterations {
                let _ = atoma_core::parser::sdf::parse_sdf_str(&content);
            }
            let total = start.elapsed();
            let avg = total / iterations;
            let mols = atoma_core::parser::sdf::parse_sdf_str(&content)?;
            let total_atoms: usize = mols.iter().map(|m| m.n_atoms()).sum();

            println!("  ┌─────────────────────────────────┐");
            println!("  │  📊 Results                     │");
            println!("  ├─────────────────────────────────┤");
            println!("  │  Molecules:{:>16}      │", mols.len());
            println!("  │  Atoms:     {:>18}      │", total_atoms);
            println!("  │  File size: {:>15.2} MB      │", size_mb);
            println!("  │  Total:     {:>18.2?}      │", total);
            println!("  │  Average:   {:>18.2?}      │", avg);
            println!("  │  Throughput:{:>15.1} atoms/ms  │",
                total_atoms as f64 / avg.as_secs_f64() / 1000.0
            );
            println!("  └─────────────────────────────────┘");
        }
        _ => {
            // PDB
            for _ in 0..warmup {
                let _ = atoma_core::parser::pdb::parse_pdb_str(&content, path);
            }
            let start = Instant::now();
            for _ in 0..iterations {
                let _ = atoma_core::parser::pdb::parse_pdb_str(&content, path);
            }
            let total = start.elapsed();
            let avg = total / iterations;
            let mol = atoma_core::parser::pdb::parse_pdb_str(&content, path)?;

            println!("  ┌─────────────────────────────────┐");
            println!("  │  📊 Results                     │");
            println!("  ├─────────────────────────────────┤");
            println!("  │  Atoms:     {:>18}      │", mol.n_atoms());
            println!("  │  File size: {:>15.2} MB      │", size_mb);
            println!("  │  Total:     {:>18.2?}      │", total);
            println!("  │  Average:   {:>18.2?}      │", avg);
            println!("  │  Throughput:{:>15.1} atoms/ms  │",
                mol.n_atoms() as f64 / avg.as_secs_f64() / 1000.0
            );
            println!("  └─────────────────────────────────┘");
        }
    }

    Ok(())
}

fn cmd_stats(path: &str) -> anyhow::Result<()> {
    let mol = parse_pdb(path)?;
    println!("{mol:#?}");
    Ok(())
}

fn cmd_extract_sequence(path: &str) -> anyhow::Result<()> {
    let fmt = detect(path);
    let mol = match fmt {
        FileFormat::Sdf | FileFormat::Mol => {
            let mols = parse_sdf(path)?;
            mols.into_iter().next().ok_or_else(|| anyhow::anyhow!("no molecules found"))?
        }
        FileFormat::Xyz => parse_xyz(path)?,
        FileFormat::MmCif => parse_mmcif(path)?,
        _ => parse_pdb(path)?,
    };

    let seqs = analysis::extract_sequences(&mol);
    if seqs.is_empty() {
        anyhow::bail!("no amino acid sequences found in this file");
    }

    let name = mol.name.as_deref().unwrap_or(
        std::path::Path::new(path).file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
    );

    print!("{}", analysis::to_fasta(&seqs, name));
    Ok(())
}

fn cmd_convert(path: &str, target: &str, output: Option<&str>) -> anyhow::Result<()> {
    let fmt = detect(path);
    let mol = match fmt {
        FileFormat::Sdf | FileFormat::Mol => {
            let mols = parse_sdf(path)?;
            mols.into_iter().next().ok_or_else(|| anyhow::anyhow!("no molecules found"))?
        }
        FileFormat::Xyz => parse_xyz(path)?,
        FileFormat::MmCif => parse_mmcif(path)?,
        _ => parse_pdb(path)?,
    };

    let result = analysis::convert_format(&mol, target)
        .ok_or_else(|| anyhow::anyhow!("unsupported target format: {target}. Supported: pdb, xyz, sdf"))?;

    match output {
        Some(out_path) => {
            std::fs::write(out_path, &result)?;
            println!("✅ Converted to {target}: {out_path} ({} atoms)", mol.n_atoms());
        }
        None => {
            print!("{result}");
        }
    }
    Ok(())
}

fn cmd_analyze(path: &str) -> anyhow::Result<()> {
    let mol = parse_pdb(path)?;
    print!("{}", analysis::analyze(&mol));
    Ok(())
}

fn truncate(s: &str, max: usize) -> &str {
    if s.len() > max {
        &s[..max]
    } else {
        s
    }
}
