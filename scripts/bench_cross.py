#!/usr/bin/env python3
"""
Cross-benchmark: molio (Rust) vs Biopython vs RDKit vs Open Babel.
Tests parsing speed and accuracy on the same PDB file.
"""
import time, sys, os, subprocess, re

PDB_FILE = os.path.join(os.path.dirname(__file__), "test_data", "pdb", "1crn.pdb")
MOLIO_BIN = os.path.join(os.path.dirname(__file__), "target", "release", "molio")
ITERATIONS = 1000
WARMUP = 20

def parse_time_to_ms(val):
    m = re.search(r'([\d.]+)\s*(µs|ms|ns|s)', str(val))
    if not m: return 0.0
    num, unit = float(m.group(1)), m.group(2)
    return {"µs": num/1000, "ms": num, "s": num*1000, "ns": num/1e6}.get(unit, 0)

def bench(name, iters, fn):
    for _ in range(WARMUP): fn()
    times, atoms = [], None
    for i in range(iters):
        t0 = time.perf_counter()
        r = fn()
        times.append((time.perf_counter() - t0) * 1000)
        if i == 0: atoms = r
    return sum(times)/len(times), min(times), atoms, len(times)

# ─── molio ────────────────────────────────────────────────────
def run_molio():
    r = subprocess.run([MOLIO_BIN, "bench", PDB_FILE, "--iterations", str(ITERATIONS)],
                       capture_output=True, text=True, timeout=30)
    atoms, avg_s, tp = None, None, None
    for line in r.stdout.split('\n'):
        if 'Atoms:' in line:
            p = line.split('Atoms:')[1].strip().split()
            if p: atoms = int(p[0])
        if 'Average:' in line: avg_s = line.split('Average:')[1].strip()
        if 'Throughput:' in line: tp = line.split('Throughput:')[1].strip().split()[0]
    return atoms, avg_s, tp

# ─── Biopython ────────────────────────────────────────────────
def run_biopython():
    from Bio.PDB.PDBParser import PDBParser
    def p():
        s = PDBParser(QUIET=True).get_structure('x', PDB_FILE)
        return sum(1 for _ in s.get_atoms())
    return bench("Bio", min(ITERATIONS, 100), p)

# ─── RDKit ────────────────────────────────────────────────────
def run_rdkit():
    from rdkit import Chem
    def p():
        m = Chem.MolFromPDBFile(PDB_FILE, removeHs=False, sanitize=False)
        return m.GetNumAtoms() if m else 0
    return bench("RDKit", min(ITERATIONS, 100), p)

# ─── Open Babel ───────────────────────────────────────────────
def run_openbabel():
    from openbabel import pybel
    def p():
        m = next(pybel.readfile("pdb", PDB_FILE))
        return len(m.atoms)
    return bench("OBabel", min(ITERATIONS, 100), p)

# ─── Raw I/O ──────────────────────────────────────────────────
def run_raw():
    def p():
        with open(PDB_FILE, 'rb') as f: return len(f.read())
    return bench("RawIO", ITERATIONS, p)

# ─── Main ─────────────────────────────────────────────────────
def main():
    print("=" * 70)
    print("  ⚛️  molio — Cross-Tool Benchmark")
    print("=" * 70)
    print(f"  File: {PDB_FILE}")
    print(f"  Size: {os.path.getsize(PDB_FILE):,} bytes | Iterations: {ITERATIONS}")
    print()

    results = {}

    print("  Running...")
    print("  [molio]      Rust native  ", end="", flush=True)
    a, s, tp = run_molio()
    ms = parse_time_to_ms(s)
    results['⚡ molio (Rust)'] = (a, ms)
    print(f"→ {ms*1000:.0f}µs | {tp} atoms/ms")

    print("  [Biopython]  PDBParser    ", end="", flush=True)
    avg, best, a, n = run_biopython()
    results['Biopython'] = (a, avg)
    print(f"→ {avg:.3f}ms | best {best:.3f}ms | {a} atoms")

    print("  [RDKit]      MolFromPDB   ", end="", flush=True)
    avg, best, a, n = run_rdkit()
    results['RDKit'] = (a, avg)
    print(f"→ {avg:.3f}ms | best {best:.3f}ms | {a} atoms")

    print("  [OpenBabel]  pybel        ", end="", flush=True)
    avg, best, a, n = run_openbabel()
    results['Open Babel'] = (a, avg)
    print(f"→ {avg:.3f}ms | best {best:.3f}ms | {a} atoms")

    print("  [Raw I/O]    disk read    ", end="", flush=True)
    avg, best, a, n = run_raw()
    results['💾 Raw disk I/O'] = (a, avg)
    print(f"→ {avg:.3f}ms | {a} bytes")

    # ─── Table ──────────────────────────────────────────────────
    molio_ms = results['⚡ molio (Rust)'][1]
    raw_ms = results['💾 Raw disk I/O'][1]
    molio_atoms = results['⚡ molio (Rust)'][0]

    print("\n" + "=" * 70)
    print("  📊 RESULTS")
    print("=" * 70)
    print(f"  {'Tool':<24} {'Atoms':>7} {'Avg Time':>12} {'vs molio':>14}")
    print(f"  {'─'*24} {'─'*7} {'─'*12} {'─'*14}")
    for name, (atoms, avg_ms) in results.items():
        a_str = str(atoms) if atoms else "—"
        if molio_ms > 0:
            ratio = avg_ms / molio_ms
            vs = f"{ratio:.0f}x slower" if ratio > 1.05 else "— BASELINE —"
        else:
            vs = "—"
        print(f"  {name:<24} {a_str:>7} {avg_ms:>9.4f}ms {vs:>14}")

    # ─── Speedups ───────────────────────────────────────────────
    print("\n" + "─" * 70)
    if molio_ms > 0:
        for name, (atoms, avg_ms) in results.items():
            if "molio" in name or "I/O" in name:
                continue
            ratio = avg_ms / molio_ms
            match = "✅ atom count match" if atoms == molio_atoms else f"⚠️ atoms: {molio_atoms} vs {atoms}"
            print(f"  🔥 molio is {ratio:.0f}x faster than {name:<12} ({match})")

    print(f"\n  📝 Parsing overhead vs raw I/O:")
    print(f"     molio:     {molio_ms/raw_ms:.1f}x  (near-I/O speed)")
    for name, (atoms, avg_ms) in results.items():
        if "molio" in name or "I/O" in name:
            continue
        print(f"     {name:<10}: {avg_ms/raw_ms:.0f}x  ({avg_ms/raw_ms:.1f}x slower than disk)")

    # ─── molio view ─────────────────────────────────────────────
    print("\n" + "─" * 70)
    r = subprocess.run([MOLIO_BIN, "view", PDB_FILE], capture_output=True, text=True)
    print(r.stdout)

if __name__ == "__main__":
    main()
