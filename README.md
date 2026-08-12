# ⚛️ atoma

**Universal Molecular File Parser & Analyzer**

Drop any molecular file → auto-detect format → parse in microseconds → instant analysis.

*Rust. WASM. Desktop. CLI. No Python. No install.*

---

## 🚀 Run it now — no install

**Web (any device, just a browser):**
```
👉 https://Dishant707.github.io/molio/web.html
```
Drop a `.pdb` / `.sdf` / `.xyz` / `.cif` file → instant analysis in your browser. 126 KB WASM engine, zero server, zero cost.

**CLI (macOS / Linux / Windows):**
```bash
cargo install atoma
atoma view 1crn.pdb
```

**Desktop GUI:**
Download `atoma-app` from [Releases](https://github.com/Dishant707/molio/releases) — 12 MB, drag-drop, no dependencies.

---

## 📊 Benchmarks vs Industry Standards

### PDB parsing (25K atoms, 100 iterations)

| Tool | Time | atoma vs |
|------|------|----------|
| **atoma** | **7.2 ms** | 1.0× |
| RDKit | 35 ms | 4.9× |
| Biopython | 90 ms | 12.5× |
| OpenBabel | 100 ms | 13.9× |

### Real-world PDB structures

| Structure | Atoms | atoma | Biopython | Speedup |
|-----------|-------|-------|-----------|---------|
| 3TAN | 11,054 | 6.3ms | 32.7ms | **5.2×** |
| 2DHB | 2,289 | 3.0ms | 7.6ms | **2.5×** |

---

## 🧬 What it does

| | |
|---|---|
| **4 formats** | PDB, mmCIF, SDF/MOL, XYZ — auto-detected |
| **Analysis** | Bonds, steric clashes, Ramachandran φ/ψ, secondary structure, Shannon entropy |
| **Sequence** | Extract amino acid sequences, export FASTA |
| **Conversion** | PDB ↔ XYZ ↔ SDF in microseconds |
| **Validation** | 28 tests, 5000 fuzz vs RDKit, Wilson CI ≥ 99.92% |

---

## 🛠️ Why atoma over existing tools?

| | atoma | Biopython | RDKit | OpenBabel |
|---|---|---|---|---|
| CLI binary | **918 KB** | 50+ MB | 200+ MB | 100+ MB |
| Desktop app | **12 MB** | ❌ | ❌ | ❌ |
| Web (WASM) | **126 KB** | ❌ | ❌ | ❌ |
| Startup | Instant | ~1s | ~2s | ~3s |
| Dependencies | **None** | Python, NumPy | Conda/Pip | apt/brew |
| Auto-detect | ✅ | ❌ | ❌ | ❌ |

---

## 📦 Install

```bash
# CLI (from crates.io)
cargo install atoma

# From source
git clone https://github.com/Dishant707/molio
cd molio && cargo build --release
./target/release/atoma --help
```

## 📖 Usage

```bash
atoma view structure.pdb          # Quick view
atoma analyze structure.pdb       # Full analysis
atoma bench structure.pdb         # Benchmark speed
atoma extract-sequence 1crn.pdb   # Extract sequence
atoma convert 1crn.pdb --to xyz   # Convert formats
```

## 🔬 Validation

- **28 unit tests** — all passing
- **5,000 differential fuzz** vs RDKit — 100% match
- **Wilson CI proof** — P(correct) ≥ 99.92%
- **Cross-tool benchmarks** — verified against Biopython, RDKit, OpenBabel

## 📄 License

MIT — use anywhere, any project.

---

*Drop a file. Get answers. No Python. No conda. No install.*
