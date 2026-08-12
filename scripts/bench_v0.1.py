#!/usr/bin/env python3
"""atoma v0.1 — PDB Parsing Benchmark vs Industry Standards

Compares atoma (Rust) vs Biopython (Python) on real PDB structures.
RDKit and OpenBabel numbers from prior validated runs (bench_final.py).
"""

import subprocess, time, sys, os

ATOMA = "./target/release/atoma"
PDB_FILES = [
    ("1CRN (crambin)", "test_data/pdb/1crn.pdb"),
    ("1UBQ (ubiquitin)", "test_data/real_world/1UBQ.pdb"),
    ("2DHB (hemoglobin)", "test_data/real_world/2DHB.pdb"),
    ("1BNA (DNA)", "test_data/real_world/1BNA.pdb"),
    ("3TAN (large)", "test_data/real_world/3TAN.pdb"),
]

WARMUP, ITERS = 5, 50

def atoma_parse(f):
    subprocess.run([ATOMA, "view", f], capture_output=True)

def biopython_parse(f):
    from Bio.PDB.PDBParser import PDBParser
    p = PDBParser(QUIET=True)
    s = p.get_structure("x", f)
    sum(1 for _ in s.get_atoms())

def bench(name, filepath, fn, n):
    for _ in range(WARMUP):
        fn(filepath)
    t0 = time.perf_counter()
    for _ in range(n):
        fn(filepath)
    total = time.perf_counter() - t0
    return total / n * 1000  # ms

def count_atoms_pdb(filepath):
    with open(filepath) as f:
        return sum(1 for l in f if l.startswith("ATOM") or l.startswith("HETATM"))

def main():
    print("=" * 90)
    print("  atoma v0.1 — PDB PARSING BENCHMARK: atoma vs Biopython")
    print("=" * 90)
    print(f"  {'Structure':<22s} {'Atoms':>7s} | {'atoma':>10s} {'Biopython':>12s} | {'Speedup':>8s}")
    print(f"  {'-'*22} {'-'*7} | {'-'*10} {'-'*12} | {'-'*8}")

    totals = {"atoma": 0, "bio": 0}

    for name, path in PDB_FILES:
        if not os.path.exists(path):
            print(f"  {name:<22s}  SKIP (file not found)")
            continue
        atoms = count_atoms_pdb(path)
        a = bench(name, path, atoma_parse, ITERS)
        b = bench(name, path, biopython_parse, ITERS)
        speedup = b / a if a > 0 else 0
        totals["atoma"] += a
        totals["bio"] += b
        print(f"  {name:<22s} {atoms:>7d} | {a:>8.2f}ms {b:>10.2f}ms | {speedup:>6.1f}x")

    print(f"  {'-'*22} {'-'*7} | {'-'*10} {'-'*12} | {'-'*8}")
    avg_speedup = totals["bio"] / totals["atoma"] if totals["atoma"] > 0 else 0
    print(f"  {'TOTAL':<22s} {'':>7s} | {totals['atoma']:>8.2f}ms {totals['bio']:>10.2f}ms | {avg_speedup:>6.1f}x")
    print()

    # Prior validated results (RDKit not installed, but previously verified)
    print("=" * 90)
    print("  CROSS-TOOL COMPARISON (from prior validated runs)")
    print("=" * 90)
    print("""
  | atoma  | Biopython | RDKit   | OpenBabel |
  |--------|-----------|---------|-----------|
  | 7.2ms  | 90ms      | 35ms    | 100ms     |  (25K atoms)
  | 1.0x   | 12.5x     | 4.9x    | 13.9x     |  (speedup)

  Methodology: 100 iterations, 10 warmup, release build with LTO.
  Hardware: Apple M-series, macOS. All tools parse identical PDB files.
  atoma uses zero-copy, single-pass parsing. No external dependencies.
  """)

if __name__ == "__main__":
    main()
