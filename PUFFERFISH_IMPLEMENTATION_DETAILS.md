# Pufferfish Implementation: Technical Deep-Dive & Performance Plan

This document provides a detailed specification for implementing the Pufferfish genome index (Almodaresi et al., 2018) and the PuffAligner (2020) extensions. It focuses on the memory-hierarchy-optimized structures and provides a concrete plan for performance measurement.

---

## 1. Data Structures: Architectural Details

The Pufferfish index is designed to minimize the number of random memory accesses by combining a Minimal Perfect Hash Function (MPHF) with a dense position table and a packed sequence array.

### A. Compacted De Bruijn Graph (cdBG) Storage
The cdBG consists of unitigs. In the "Flat" representation:
- **Sequence Array (S)**: All unitigs are concatenated into a single bit-packed array. For 2-bit packing (A=00, C=01, G=10, T=11), 4 bases fit in 1 byte.
- **Unitig Boundaries**: Instead of a separate array of pointers, we use the **Blocked Layout** (as specified in `FLAT_INDEX_PLAN.md`). Each 64-byte block (matching the L1 cache line) is self-describing.

### B. MPHF (The Entry Point)
- **Role**: Maps each unique canonical k-mer to a unique integer in the range `[0, N-1]`.
- **Implementation**: Uses a function like **BBHash** or **PTHash**. It must be "minimal" (no holes in the output range) and "perfect" (no collisions).
- **Cache Impact**: The MPHF itself is a data structure (e.g., bit-vectors and rank/select structures) that must be consulted. It usually requires 3–4 bits per k-mer.

### C. The POS Table (The Mapping)
- **Structure**: A dense array of 64-bit integers where `index = MPHF(kmer)`.
- **Packing Scheme**:
  - `[0-55]` bits: **Block Offset** (Absolute byte offset in the Flat Array).
  - `[56-63]` bits: **Local K-mer Offset** (Index of the k-mer within that specific 64-byte block).
- **Optimization**: This table is accessed once per seed k-mer. It is the most frequent source of LLC (Last Level Cache) misses.

---

## 2. The Extension Algorithm (PuffAligner)

Once a seed k-mer is located via the POS table, subsequent bases are aligned using the **Extension** method, which avoids further hash lookups.

### Algorithm: Seed-and-Extend
1.  **Seed**: Look up the first k-mer of the read in the MPHF/POS table.
2.  **Locate**: Fetch the 64-byte block from the Flat Array.
3.  **Local Match**:
    - Compare the read's next base with the next base in the block's packed DNA.
    - Since 4 bases are packed per byte, use bit-shifts to extract the next 2 bits.
4.  **Boundary Jump**:
    - If the end of the unitig is reached, consult the **Edge Table (ETAB)** within the *same* 64-byte block.
    - Find the edge corresponding to the next base in the read.
    - Update the current pointer to the `to_block_offset` stored in the ETAB row.
5.  **Termination**: If no edge matches or a mismatch occurs in the DNA, stop the alignment.

---

## 3. Performance Measurement Plan

To satisfy the requirements of the project, we must quantify the efficiency of this layout compared to the Flat List approach.

### A. Measuring Cache Misses
- **Tool**: `perf` on Linux or `xctrace` (Instruments) on macOS.
- **Metrics**:
  - `L1-dcache-load-misses`: Measures spatial locality. High misses suggest poor data grouping.
  - `LLC-load-misses` (Last Level Cache): Measures temporal locality and memory bandwidth saturation. This is the "Gold Standard" for indexing performance.
- **Procedure**:
  ```bash
  # Linux example
  perf stat -e L1-dcache-load-misses,LLC-load-misses ./bench queries.fasta
  ```

### B. Runtime vs. Flat List
- **Baseline (Flat List)**: A `std::unordered_map<string, vector<Occurrence>>`.
- **Comparison Points**:
  - **Latency**: Time per single k-mer query (nanoseconds).
  - **Throughput**: Millions of queries per second (M ops/s).
  - **Effect of Query Count**: Sweep from 10k to 10M queries to see when the index exceeds the cache size.

### C. Memory Usage
- **Metric**: Bits per k-mer (bpk).
- **Calculation**:
  - `Total Memory = (Size of MPHF + Size of POS Table + Size of Flat Array)`.
  - `BPK = (Total Memory * 8) / Total unique k-mers`.
- **Target**: Pufferfish typically achieves < 15-20 bpk, whereas Flat List maps often use > 100 bpk.

---

## 4. Implementation Checklist: Cache-Awareness

- [ ] **Alignment**: Ensure the Flat Array is aligned to 64-byte boundaries (`std::aligned_alloc`).
- [ ] **Prefetching**: (Optional) Experiment with `__builtin_prefetch` for the next k-mer's POS entry while processing the current one.
- [ ] **Bit-Packing**: Use `uint64_t` shifts and masks instead of `std::bitset` to ensure the compiler generates efficient SIMD-friendly code.
- [ ] **Blocked Layout**: All metadata (length, edges, occurrences) for a unitig must be reachable within one or two cache-line fetches from the start of the unitig's block.

---

## 5. Benchmarking Deliverables
At the end of the implementation, the following table must be filled in `bench-report.md`:

| Metric | Flat List | Pufferfish |
| :--- | :--- | :--- |
| **Throughput (M ops/s)** | | |
| **L1 Misses / Query** | | |
| **LLC Misses / Query** | | |
| **Memory (Bits/K-mer)** | | |
| **Build Time (s)** | | |
