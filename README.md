# Flat Index — Cache-Efficient Genome Indexing

A compacted de Bruijn graph index that packs DNA, occurrence data (UTAB), and edge metadata (ETAB) into 64-byte cache-line-aligned blocks for efficient k-mer lookup and read alignment.

## Build

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j$(sysctl -n hw.logicalcpu)
cd ..
```

## Usage

See [WORKFLOW.md](WORKFLOW.md) for full pipeline instructions (two paths: project's own builder or Cuttlefish if you can get it to work).

Quick start with the project builder:

```bash
# 1. Build the compacted de Bruijn graph
./build/build_unitigs genome.fasta data/genome 31

# 2. Build the occurrence table
./build/build_occ genome.fasta data/genome 31

# 3. Build the flat index
./build/flatindex_build --tsv data/genome --k 31 --out genome.flat

# 4. Query
./build/flatindex_query genome.flat --kmer ACGTACGTACGTACGTACGTACGTACGTACG
./build/flatindex_query genome.flat --fastq reads.fq
```

## Benchmark

Compares flat index against Pufferfish on alignment throughput. Requires Pufferfish (see [WORKFLOW.md](WORKFLOW.md) for install instructions).

```bash
bash scripts/benchmark_pufferfish.sh testfasta/ecoli_k12_sub_mg1655.fasta
```

