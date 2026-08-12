#!./.venv/bin/python3
"""
Differential fuzzing: generates random molecules, parses with BOTH
molio and RDKit, compares results. Any mismatch is a potential bug.
"""
import subprocess, sys, os, random, math, tempfile

MOLIO_BIN = os.path.join(os.path.dirname(__file__), "target", "release", "molio")
N_TESTS = 5000
ELEMENTS = ["C", "N", "O", "S", "P", "F", "Cl", "Br"]

def gen_random_xyz(n_atoms):
    """Generate a random XYZ file."""
    lines = [f"{n_atoms}", f"Fuzz molecule {random.randint(1, 99999)}"]
    for _ in range(n_atoms):
        elem = random.choice(ELEMENTS)
        x = random.uniform(-50, 50)
        y = random.uniform(-50, 50)
        z = random.uniform(-50, 50)
        lines.append(f"{elem}  {x:.6f}  {y:.6f}  {z:.6f}")
    return "\n".join(lines) + "\n"

def gen_random_sdf(n_atoms):
    """Generate a random SDF file."""
    n_bonds = max(0, n_atoms - 1)
    lines = ["FuzzMol", "  generated", ""]
    lines.append(f"{n_atoms:>3}{n_bonds:>3}  0  0  0  0  0  0  0  0999 V2000")
    for _ in range(n_atoms):
        elem = random.choice(ELEMENTS)
        x = random.uniform(-30, 30)
        y = random.uniform(-30, 30)
        z = random.uniform(-30, 30)
        lines.append(f"{x:>10.4f}{y:>10.4f}{z:>10.4f} {elem:<3} 0  0  0  0  0  0  0  0  0  0  0  0")
    for i in range(n_bonds):
        a1 = i + 1
        a2 = min(i + 2, n_atoms)
        if a1 != a2:
            lines.append(f"{a1:>3}{a2:>3}  1  0  0  0  0")
    lines.append("M  END\n$$$$")
    return "\n".join(lines) + "\n"

def molio_parse_sdf(content):
    """Parse SDF with molio via temp file, return (n_mols, total_atoms)."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.sdf', delete=False) as f:
        f.write(content)
        tmp = f.name
    try:
        r = subprocess.run([MOLIO_BIN, "view", tmp], capture_output=True, text=True, timeout=5)
        # Parse output: "Molecules: X" and "Atom count"
        n_mols = 0
        total = 0
        for line in r.stdout.split('\n'):
            if 'Atoms:' in line:
                p = line.split('Atoms:')[1].strip().split()
                if p:
                    total = max(total, int(p[0]))
        if 'Molecule 1/' in r.stdout:
            # Multi-molecule output
            for line in r.stdout.split('\n'):
                if 'Total:' in line:
                    parts = line.split()
                    for i, p in enumerate(parts):
                        if p == 'molecules':
                            n_mols = int(parts[i-1])
        else:
            n_mols = 1
        return n_mols, total
    except:
        return 0, 0
    finally:
        os.unlink(tmp)

def rdkit_parse_sdf(content):
    """Parse SDF with RDKit from string, return (n_mols, total_atoms)."""
    from rdkit import Chem
    import io
    suppl = Chem.SDMolSupplier()
    suppl.SetData(content, sanitize=False, removeHs=False)
    mols = [m for m in suppl if m is not None]
    n_atoms = sum(m.GetNumAtoms() for m in mols)
    return len(mols), n_atoms

def main():
    print("=" * 60)
    print("  Differential Fuzzing: molio vs RDKit")
    print("=" * 60)
    print(f"  Tests: {N_TESTS} random molecules")
    print()

    mismatches = 0
    rdkit_fails = 0
    molio_fails = 0

    for i in range(N_TESTS):
        n_atoms = random.randint(1, 20)
        sdf = gen_random_sdf(n_atoms)

        # Parse with both
        try:
            rd_n, rd_atoms = rdkit_parse_sdf(sdf)
        except Exception as e:
            rdkit_fails += 1
            continue

        try:
            ml_n, ml_atoms = molio_parse_sdf(sdf)
        except:
            molio_fails += 1
            continue

        # Compare
        if rd_atoms != ml_atoms:
            mismatches += 1
            if mismatches <= 5:  # Show first 5 mismatches only
                print(f"  ⚠️  MISMATCH #{mismatches}:")
                print(f"     molio: {ml_atoms} atoms | RDKit: {rd_atoms} atoms")
                print(f"     SDF snippet: {sdf[:100]}...")

        if (i + 1) % 1000 == 0:
            print(f"  ... {i+1}/{N_TESTS} tested, {mismatches} mismatches so far")

    print()
    print("=" * 60)
    print("  RESULTS")
    print("=" * 60)
    print(f"  Total tests:      {N_TESTS}")
    print(f"  Mismatches:       {mismatches} ({100*mismatches/N_TESTS:.2f}%)")
    print(f"  RDKit failures:   {rdkit_fails}")
    print(f"  molio failures:   {molio_fails}")
    print(f"  Agreement rate:   {100*(1 - mismatches/(N_TESTS-rdkit_fails-molio_fails)):.2f}%")
    print()

    if mismatches == 0:
        print("  ✅ molio and RDKit agree on all random molecules!")
    else:
        print(f"  ⚠️  {mismatches} mismatches found — investigate above")

if __name__ == "__main__":
    main()
