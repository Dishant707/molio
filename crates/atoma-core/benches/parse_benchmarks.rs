//! Criterion benchmarks for molio-core parsers.
//!
//! Run with: `cargo bench -p molio-core`

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use atoma_core::parser::pdb::parse_pdb_str;

const PDB_1CRN: &str = include_str!("../../test_data/pdb/1crn.pdb");

fn bench_parse_pdb(c: &mut Criterion) {
    c.bench_function("parse_pdb_1crn", |b| {
        b.iter(|| parse_pdb_str(black_box(PDB_1CRN), black_box("1crn.pdb")))
    });
}

fn bench_parse_pdb_large(c: &mut Criterion) {
    // Parse 1CRN 100 times to simulate a larger structure
    let large = PDB_1CRN.repeat(100);
    c.bench_function("parse_pdb_1crn_x100", |b| {
        b.iter(|| parse_pdb_str(black_box(&large), black_box("large.pdb")))
    });
}

criterion_group!(benches, bench_parse_pdb, bench_parse_pdb_large);
criterion_main!(benches);
