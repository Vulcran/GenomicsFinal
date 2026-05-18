# Flat Index — End-to-End Workflows

Two paths to build and query the flat index from real genomic data.

---

## Path A — Project Pipeline (build_unitigs + build_occ)

Uses the project's own cdBG builder. Good for development and correctness testing.



mkdir build
cd build
cmake ..
cmake --build . --config Release
cd ..


```
genome.fasta
     │
     ▼
bin/build_unitigs   →   <prefix>.unitigs.tsv
                    →   <prefix>.edges.tsv
     │
     ▼
bin/build_occ       →   <prefix>.occ.tsv
                    →   <prefix>.refs.tsv
     │
     ▼
flatindex_build     →   genome.flat
     │
     ▼
flatindex_query     →   hits
```

### Step-by-step

```bash
# 1. Build the compacted de Bruijn graph
./Project/bin/build_unitigs genome.fasta data/genome 31
# writes: data/genome.unitigs.tsv  data/genome.edges.tsv

# 2. Build the occurrence table
./Project/bin/build_occ genome.fasta data/genome 31
# writes: data/genome.occ.tsv  data/genome.refs.tsv

# 3. Build the flat index
./build/flatindex_build --tsv data/genome --k 31 --out genome.flat

# 4. Query
./build/flatindex_query genome.flat --kmer ACGTACGTACGTACGTACGTACGTACGTACG
./build/flatindex_query genome.flat --fastq reads.fq

# 5. Benchmark vs naive baseline
./build/bench --tsv data/genome --k 31
```

---

## Path B — Cuttlefish

Uses [Cuttlefish 2](https://github.com/COMBINE-lab/cuttlefish) as the cdBG builder.
Cuttlefish is faster and scales to full mammalian genomes; the occurrence table is
still computed by `build_occ` (which only needs the unitig sequences, not how they
were built).

```
genome.fasta
     │
     ▼
cuttlefish build    →   output.gfa   (segments + links)
     │
     ├── awk/python GFA→TSV conversion
     │
     ▼
bin/build_occ       →   output.occ.tsv
                    →   output.refs.tsv
     │
     ▼
flatindex_build     →   genome.flat
  (--gfa + --occ)
     │
     ▼
flatindex_query     →   hits
```

### Step-by-step

**Install Cuttlefish** (pick one):
```bash
conda install -c bioconda cuttlefish   # via Bioconda

or if ths doesnt work

CONDA_SUBDIR=osx-64 conda install -c bioconda cuttlefish
# or build from source:
git clone https://github.com/COMBINE-lab/cuttlefish && cd cuttlefish
mkdir build && cd build && cmake .. && make -j8




```

**Run the pipeline:**
```bash
# 1. Build the cdBG with Cuttlefish (GFA output includes segments + links)
# or run conda activate base before. then you dont need conda run
conda run cuttlefish build \
  -r /testfasta/200bp.fasta \
  -k 31 \
  -t 8 \
  -o output \
  -f 3 \
  -w /tmp/cf_work
# -f 3 = GFA-reduced (segments + links only)
# -w sets a temp working directory for KMC k-mer counting
# writes: output.gfa

# 2. Convert GFA segments to unitigs.tsv so build_occ can read them.
#    IDs are assigned 0-indexed in segment-appearance order — this matches
#    what flatindex_build's GFA loader does internally.
awk 'BEGIN{print "id\tlength\tsequence"; n=0}
     /^S/  {print n"\t"length($3)"\t"$3; n++}' \
  output.gfa > output.unitigs.tsv

# 3. Build the occurrence table from the reference + Cuttlefish unitigs
./Project/bin/build_occ genome.fasta output 31
# writes: output.occ.tsv  output.refs.tsv

# 4. Build the flat index (GFA for unitigs+edges, TSV for occurrences)
./build/flatindex_build \
  --gfa output.gfa \
  --occ output.occ.tsv \
  --k 31 \
  --out genome.flat

# 5. Query
./build/flatindex_query genome.flat --kmer ACGTACGTACGTACGTACGTACGTACGTACG
./build/flatindex_query genome.flat --fastq reads.fq
```

---

## Comparison

| | Path A (build_unitigs) | Path B (Cuttlefish) |
|---|---|---|
| cdBG builder | Project's own (C++20, in-repo) | Cuttlefish 2 (external, faster) |
| Genome scale | Up to ~100 MB unique k-mers | Full mammalian genomes |
| Occurrence table | `build_occ` | `build_occ` (same tool) |
| Index build | `flatindex_build --tsv` | `flatindex_build --gfa --occ` |
| Query | `flatindex_query` | `flatindex_query` (identical) |
| Extra install | None | Cuttlefish (conda or source) |

The flat index and query code are identical in both paths. Only the input source changes.

---

## Benchmark with real data

```bash
# Run the throughput + cache-miss comparison on real data
./build/bench --tsv data/genome --k 31

# On Linux, wrap with perf for hardware cache counters:
perf stat -e L1-dcache-load-misses,LLC-load-misses \
  ./build/bench --tsv data/genome --k 31

# On macOS, use xctrace:
xctrace record --template "CPU Counters" \
  --launch -- ./build/bench --tsv data/genome --k 31
```
