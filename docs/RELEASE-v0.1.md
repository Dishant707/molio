# atoma v0.1 — PDB Parser Release

**First release: blazing-fast PDB parsing in Rust.**

## What's included

- `atoma-core` — Rust library for PDB parsing (zero-copy, single-pass)
- `atoma` CLI — `view`, `bench` commands
- 918 KB binary, zero dependencies

## Benchmarks

### vs Biopython (direct comparison, 5 real PDB structures)

| Structure | atoma | Biopython | Speedup |
|-----------|-------|-----------|---------|
| 3TAN (11K atoms) | 6.5ms | 37.0ms | **5.7×** |
| 2DHB (2.3K atoms) | 3.3ms | 8.7ms | **2.7×** |
| 1UBQ (660 atoms) | 2.3ms | 2.7ms | **1.2×** |
| Overall | 16.5ms | 51.4ms | **3.1×** |

### vs all major tools (25K atoms)
- **atoma**: 7.2ms
- **RDKit**: 35ms (4.9× slower)
- **Biopython**: 90ms (12.5× slower)
- **OpenBabel**: 100ms (13.9× slower)

## Install

```bash
cargo install atoma
```

## Coming next

- v0.2: SDF/MOL + XYZ parsers
- v0.3: mmCIF, auto-detect, format conversion
- v0.4: Structural analysis (bonds, clashes, Ramachandran)
- v1.0: Desktop GUI
