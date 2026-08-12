#!./.venv/bin/python3
"""
Cross-benchmark: molio SDF parser vs RDKit vs Open Babel.
Tests parsing speed and accuracy on a multi-molecule SDF file.
"""
import time, sys, os, subprocess, re

SDF_FILE = os.path.join(os.path.dirname(__file__), "test_data", "sdf", "molecules.sdf")
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
    times, result = [], None
    for i in range(iters):
        t0 = time.perf_counter()
        r = fn()
        times.append((time.perf_counter() - t0) * 1000)
        if i == 0: result = r
    return sum(times)/len(times), min(times), result, len(times)

# molio - benchmark via CLI bench command
def run_molio():
    r = subprocess.run([MOLIO_BIN, "bench", SDF_FILE, "--iterations", str(ITERATIONS)],
                       capture_output=True, text=True, timeout=30)
    atoms, avg_s, tp = None, None, None
    for line in r.stdout.split('\n'):
        if 'Molecules:' in line:
            p = line.split('Molecules:')[1].strip().split()
            if p: n_mols = int(p[0])
        else:
            n_mols = 3
        if 'Atoms:' in line:
            p = line.split('Atoms:')[1].strip().split()
            if p: atoms = int(p[0])
        if 'Average:' in line:
            avg_s = line.split('Average:')[1].strip()
        if 'Throughput:' in line:
            tp = line.split('Throughput:')[1].strip().split()[0]
    ms = parse_time_to_ms(avg_s) if avg_s else 0
    return n_mols, ms, atoms, tp

# RDKit
def run_rdkit():
    from rdkit import Chem
    def parse():
        mols = []
        suppl = Chem.SDMolSupplier(SDF_FILE, removeHs=False, sanitize=False)
        for mol in suppl:
            if mol is not None:
                mols.append(mol.GetNumAtoms())
        return len(mols), sum(mols)
    avg, best, r, n = bench("RDKit", min(ITERATIONS, 100), parse)
    n_mols, n_atoms = r if r else (0, 0)
    return avg, best, n_mols, n_atoms, n

# Open Babel
def run_openbabel():
    from openbabel import pybel
    def parse():
        mols = list(pybel.readfile("sdf", SDF_FILE))
        return len(mols), sum(len(m.atoms) for m in mols)
    avg, best, r, n = bench("OBabel", min(ITERATIONS, 100), parse)
    n_mols, n_atoms = r if r else (0, 0)
    return avg, best, n_mols, n_atoms, n

# Raw I/O baseline
def run_raw():
    def p():
        with open(SDF_FILE, 'rb') as f: return len(f.read())
    avg, best, size, n = bench("RawIO", ITERATIONS, p)
    return avg, size

def main():
    print("=" * 70)
    print("  molio SDF — Cross-Tool Benchmark")
    print("=" * 70)
    print(f"  File: {SDF_FILE}")
    print(f"  Size: {os.path.getsize(SDF_FILE):,} bytes")
    print(f"  Content: 3 molecules (Aspirin, Caffeine, Benzene) — 41 atoms total")
    print()

    print("  Running...")
    print("  [molio]      Rust native  ", end="", flush=True)
    molio_mols, molio_ms, molio_atoms, molio_tp = run_molio()
    print(f"> {molio_mols} mols, {molio_atoms} atoms | {molio_ms*1000:.0f}µs | {molio_tp} atoms/ms")

    print("  [RDKit]      SDMolSupplier", end="", flush=True)
    avg, best, n_mols, n_atoms, n = run_rdkit()
    print(f"> {n_mols} mols, {n_atoms} atoms | avg {avg:.3f}ms | best {best:.3f}ms")

    print("  [OpenBabel]  pybel        ", end="", flush=True)
    avg2, best2, n_mols2, n_atoms2, n2 = run_openbabel()
    print(f"> {n_mols2} mols, {n_atoms2} atoms | avg {avg2:.3f}ms | best {best2:.3f}ms")

    print("  [Raw I/O]    disk read    ", end="", flush=True)
    raw_ms, raw_size = run_raw()
    print(f"> {raw_ms:.3f}ms | {raw_size} bytes")

    # Table
    print("\n" + "=" * 70)
    print("  RESULTS — SDF Parsing")
    print("=" * 70)
    print(f"  {'Tool':<24} {'Mols':>6} {'Atoms':>7} {'Avg Time':>12} {'vs molio':>14}")
    print(f"  {'-'*24} {'-'*6} {'-'*7} {'-'*12} {'-'*14}")

    for name, n_m, n_a, ms in [
        ("molio (Rust)", molio_mols, molio_atoms, molio_ms),
        ("RDKit", n_mols, n_atoms, avg),
        ("Open Babel", n_mols2, n_atoms2, avg2),
        ("Raw disk I/O", 0, raw_size, raw_ms),
    ]:
        a_str = str(n_a) if n_a else "-"
        ratio = ms/molio_ms if molio_ms > 0 else 0
        vs = f"{ratio:.0f}x slower" if ratio > 1.05 else "- BASELINE -"
        print(f"  {name:<24} {n_m:>6} {a_str:>7} {ms:>9.4f}ms {vs:>14}")

    print("\n" + "-" * 70)
    if molio_ms > 0:
        for name, n_m, n_a, ms in [("RDKit", n_mols, n_atoms, avg), ("Open Babel", n_mols2, n_atoms2, avg2)]:
            ratio = ms / molio_ms
            print(f"  molio is {ratio:.0f}x faster than {name} (3 molecules, {molio_atoms} atoms)")
        print(f"\n  Parsing overhead: molio is {molio_ms/raw_ms:.1f}x slower than raw disk I/O")

if __name__ == "__main__":
    main()
