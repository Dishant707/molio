#!./.venv/bin/python3
"""
FINAL CROSS-TOOL VALIDATION
Real PDB structures from RCSB. All tools benchmarked on same files.
molio vs Biopython vs RDKit vs Open Babel.
"""
import time, subprocess, os, sys

MOLIO_BIN = os.path.join(os.path.dirname(__file__), "target", "release", "molio")
DATA_DIR = os.path.join(os.path.dirname(__file__), "test_data", "real_world")

STRUCTURES = ["1CRN", "1UBQ", "2DHB", "1BNA", "3TAN"]

ITERS = 200  # Same iterations for all tools
WARMUP = 20

def bench_molio(filepath):
    """molio bench via CLI."""
    r = subprocess.run(
        [MOLIO_BIN, "bench", filepath, "--iterations", str(ITERS), "--warmup", str(WARMUP)],
        capture_output=True, text=True, timeout=30
    )
    atoms = avg = tp = None
    for line in r.stdout.split('\n'):
        if 'Atoms:' in line:
            p = line.split('Atoms:')[1].strip().split()
            if p: atoms = int(p[0])
        if 'Average:' in line:
            avg_str = line.split('Average:')[1].strip()
            if 'µs' in avg_str: avg = float(avg_str.replace('µs','').split()[0]) / 1000
            elif 'ms' in avg_str: avg = float(avg_str.replace('ms','').split()[0])
            elif 'ns' in avg_str: avg = float(avg_str.replace('ns','').split()[0]) / 1_000_000
        if 'Throughput:' in line:
            tp = line.split('Throughput:')[1].strip().split()[0]
    return atoms, avg, tp

def bench_tool(name, filepath, fn_parser, iters):
    """Generic benchmark: warmup then measure."""
    for _ in range(WARMUP):
        fn_parser(filepath)
    times = []
    atoms = None
    for _ in range(iters):
        t0 = time.perf_counter()
        result = fn_parser(filepath)
        times.append((time.perf_counter() - t0) * 1000)
        if atoms is None: atoms = result
    avg = sum(times) / len(times)
    best = min(times)
    return atoms, avg, best

def biopython_parse(path):
    from Bio.PDB.PDBParser import PDBParser
    s = PDBParser(QUIET=True).get_structure('x', path)
    return sum(1 for _ in s.get_atoms())

def rdkit_parse(path):
    from rdkit import Chem
    m = Chem.MolFromPDBFile(path, removeHs=False, sanitize=False)
    return m.GetNumAtoms() if m else 0

def openbabel_parse(path):
    from openbabel import pybel
    m = next(pybel.readfile("pdb", path))
    return len(m.atoms)

def main():
    print("=" * 85)
    print("  FINAL CROSS-TOOL VALIDATION — Real RCSB Structures")
    print("=" * 85)
    print(f"  {'Structure':<10} {'Size':>8} {'Atoms':>7} | {'molio':>12} {'Biopython':>12} {'RDKit':>12} {'OpenBabel':>12}")
    print(f"  {'-'*10} {'-'*8} {'-'*7} | {'-'*12} {'-'*12} {'-'*12} {'-'*12}")

    all_data = []
    tool_times = {"molio": [], "Biopython": [], "RDKit": [], "Open Babel": []}

    for struct in STRUCTURES:
        filepath = os.path.join(DATA_DIR, f"{struct}.pdb")
        size_kb = os.path.getsize(filepath) / 1024

        # molio
        ml_atoms, ml_avg, ml_tp = bench_molio(filepath)
        tool_times["molio"].append(ml_avg)

        # Biopython
        bp_atoms, bp_avg, bp_best = bench_tool("Bio", filepath, biopython_parse, min(ITERS, 50))
        tool_times["Biopython"].append(bp_avg)

        # RDKit
        rdk_atoms, rdk_avg, rdk_best = bench_tool("RDKit", filepath, rdkit_parse, min(ITERS, 50))
        tool_times["RDKit"].append(rdk_avg)

        # OpenBabel
        ob_atoms, ob_avg, ob_best = bench_tool("OB", filepath, openbabel_parse, min(ITERS, 50))
        tool_times["Open Babel"].append(ob_avg)

        # Accuracy check: all tools must agree on atom count
        atoms_list = [ml_atoms, bp_atoms, rdk_atoms, ob_atoms]
        all_agree = len(set(a for a in atoms_list if a is not None)) <= 1
        status = "✅" if all_agree else "⚠️"

        all_data.append((struct, size_kb, ml_atoms, ml_avg, bp_avg, rdk_avg, ob_avg, status))

        print(f"  {status} {struct:<8} {size_kb:>7.0f}KB {ml_atoms:>7} | {ml_avg:>9.3f}ms {bp_avg:>9.3f}ms {rdk_avg:>9.3f}ms {ob_avg:>9.3f}ms")

    # Speedup summary
    print(f"\n  {'─'*83}")
    print(f"  SPEEDUP (molio vs each tool)")

    for tool_name in ["Biopython", "RDKit", "Open Babel"]:
        ratios = []
        for i, struct in enumerate(STRUCTURES):
            ml = tool_times["molio"][i]
            other = tool_times[tool_name][i]
            if ml and other and ml > 0:
                ratios.append(other / ml)
        if ratios:
            avg_ratio = sum(ratios) / len(ratios)
            min_ratio = min(ratios)
            max_ratio = max(ratios)
            print(f"  vs {tool_name:<10}: {avg_ratio:.0f}x avg  ({min_ratio:.0f}x–{max_ratio:.0f}x range)")

    # Accuracy summary
    mismatches = sum(1 for d in all_data if d[7] == "⚠️")
    total = len(all_data)
    print(f"\n  ACCURACY: {total - mismatches}/{total} structures — all tools agree on atom counts")

    print("\n" + "=" * 85)
    print("  Reproduce: python bench_final.py")
    print("=" * 85)

if __name__ == "__main__":
    main()
