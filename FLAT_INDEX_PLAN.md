# Flat Index — Implementation Plan


USE CUTTLEFISH?


A step-by-step plan to build the cache-friendly Pufferfish-style flat index
described in the project proposal. Scoped to the **packing + query** half of the
project; the cdBG construction half (unitigs, edges) is treated as a black-box
input that we either hardcode or import from the existing `Project/` pipeline.

The plan is ordered so each step produces a runnable artifact and can be
benchmarked independently. Earlier steps have no dependency on later ones.

---

## Glossary (so the steps below can reference these names cleanly)

| Term | Meaning |
|---|---|
| **canonical k-mer** | `min(K, RC(K))` — the lexicographically smaller of a k-mer and its reverse complement. The thing the MPHF actually hashes. |
| **MPHF** | Minimal Perfect Hash Function. Maps the `N` canonical k-mers to slots `[0, N)` with no collisions. |
| **POS table** | Dense array, one 64-bit word per MPHF slot. Bit-packed: 56-bit offset into the Flat Array + 8-bit local k-mer offset within that unitig. |
| **Flat Array** | A `std::vector<std::byte>` (or equivalent) that stores all unitig blocks back-to-back. Indexed by *byte offset*, not unitig id. |
| **Block** | A region of the Flat Array (one or more 64-byte cache lines) containing one unitig's payload: header + DNA + UTAB rows + ETAB rows. |
| **UTAB row** | "This unitig appears at FASTA position X in orientation Y" record. Comes from `build_occ`'s `occ.tsv`. |
| **ETAB row** | "This unitig connects to neighbor unitig U via base X" record. Comes from `build_unitigs`'s `edges.tsv`, augmented with the extension base at pack time. |

---

## Step 0 — Decide the block layout

**Goal**: lock down the byte layout of one block before writing any code, because every other step depends on this contract.

A 64-byte block has roughly this shape:

```text
offset  size  field
──────  ────  ──────────────────────────────────────────
 0       2    unitig_id            (uint16, debug only)
 2       1    dna_len_kmers        (uint8, # k-mers in this block)
 3       1    n_utab               (uint8, # UTAB rows)
 4       1    n_etab               (uint8, # ETAB rows)
 5       1    flags                (uint8, e.g. has_continuation)
 6       2    next_block_offset/64 (uint16, for spillover unitigs)
 8       D    DNA, 2-bit packed    (D bytes, see sizing below)
 8+D    8·U   UTAB rows            (8 bytes/row × U rows)
 ...    4·E   ETAB rows            (4 bytes/row × E rows)
 ...   pad    zero padding to 64 B
```

**UTAB row (8 bytes)** for one occurrence:

| field | bits | notes |
|---|---|---|
| ref_id    | 16 | up to 65k references |
| ref_pos   | 32 | up to 4 GB per reference |
| entry_off | 8  | offset into unitig of run start |
| walk_len  | 7  | up to 127 consecutive k-mers (fine for blocks ≤ 256 k-mers) |
| orient    | 1  | + or − |

**ETAB row (4 bytes)** for one outgoing edge:

| field | bits | notes |
|---|---|---|
| to_block_offset / 64 | 24 | up to 16M blocks (≈1 GB of blocks) |
| ext_base             | 2  | A/C/G/T |
| from_orient          | 1  | which end of *this* unitig |
| to_orient            | 1  | which end of *that* unitig |
| reserved             | 4  | padding |

**Sizing the DNA region**: with 2-bit packing, `D` bytes hold `4·D` bases =
`4·D − k + 1` k-mers. To leave ~16 bytes for ≤2 UTAB and ≤4 ETAB rows in the
common case, target `D ≈ 40 bytes` → up to 160 bp → up to 130 k-mers per
block (k=31). That fits the 8-bit local-offset constraint.

**Spillover for long unitigs**: a unitig longer than one block's DNA capacity
chains into a continuation block via `next_block_offset/64`. The continuation
block has the same header but typically zero UTAB/ETAB rows (the DNA carries
into the next block, but the metadata stayed in the head block). An MPHF slot
points at *whichever block contains its k-mer*, so the 8-bit local offset
indexes into that block's DNA, not the unitig as a whole.

**Action**: Write this layout into a single header file `flat_block.hpp` with
`constexpr` offsets, sizes, and pack/unpack helpers. Every other step should
include this header and never hardcode magic numbers.

---

## Step 1 — Hardcoded fixture and core types

**Goal**: bring up the rest of the pipeline against a synthetic example so we can
iterate on the *index* code without depending on the cdBG builder.

**Action**:

1. Create `data_model.hpp` defining plain-C++ structs for the inputs:

   ```cpp
   struct InUnitig    { uint32_t id; std::string seq; };
   struct InOccurrence{ uint32_t unitig_id; uint16_t ref_id;
                        uint32_t ref_pos; uint8_t entry_off;
                        uint8_t walk_len; char orient; };
   struct InEdge      { uint32_t from_id; char from_orient;
                        uint32_t to_id;   char to_orient; uint8_t overlap; };
   ```

2. Add a `fixtures/tiny.cpp` (or in-`main` helper) that returns hardcoded
   vectors of those three types describing one obvious example:
   - 3 unitigs, 1 reference, k = 7 (so it's hand-checkable)
   - At least one orient='-' occurrence, at least one split edge.

3. Confirm round-trip: print the fixture, eyeball it, sanity-check by hand.

**Deliverable**: `bin/flatindex_demo` that links the fixture and prints "loaded N
unitigs, M occurrences, E edges" — nothing else yet.

---

## Step 2 — Canonical k-mer toolkit

**Goal**: shared utilities used by every subsequent step. No flat array yet.

**Action**: Add `kmer.hpp` with:

- `uint64_t encode(const std::string& kmer)` — 2-bit pack into a uint64.
- `uint64_t reverse_complement(uint64_t kmer, int k)` — bit-twiddle, identical
  to `Project/scripts/build_unitigs.cpp`.
- `uint64_t canonical(uint64_t kmer, int k)` — `std::min(kmer, RC(kmer))`.
- `bool is_canonical(uint64_t kmer, int k)` — `kmer == canonical(...)`.
- `std::string decode(uint64_t kmer, int k)` — for debugging/printing.

Unit-test by canonicalizing known pairs:
`canonical("AGTAA") == canonical("TTACT") == "AGTAA"` (because AGTAA is
lexicographically smaller as a 2-bit integer).

**Deliverable**: a `tests/test_kmer.cpp` that runs ~10 hand-computed assertions.

---

## Step 3 — DNA 2-bit packer for the block

**Goal**: pack a unitig's bases into bytes the same way the block layout
expects, with offset-aware reads (so the query path can extract a k-mer at any
local offset).

**Action**: Add `dna_pack.hpp` with:

- `void pack_dna(std::string_view bases, std::byte* dst, std::size_t dst_bits)`
- `uint64_t read_kmer(const std::byte* src, std::size_t bit_off, int k)`
  — extract a 2-bit-packed k-mer starting at `bit_off` bits into `src`. Use
  shifts and a `__builtin_bswap64` if you want to be fast; correctness first.

Unit-test: pack `"ACGT"` → `0b00011011`, read it back at every offset.

**Deliverable**: `tests/test_dna_pack.cpp`.

---

## Step 4 — Block builder

**Goal**: take one unitig + its UTAB rows + its ETAB rows, return a sequence of
filled 64-byte blocks (one or more, depending on length).

**Action**: Add `block_builder.hpp` / `.cpp` with:

```cpp
struct PackedBlocks {
    std::vector<std::byte> bytes;            // contiguous, multiple of 64
    uint32_t  first_block_byte_offset;       // where in the global flat array
    uint8_t   kmers_per_block_dna;           // for the local-offset arithmetic
};

PackedBlocks pack_unitig(const InUnitig& u,
                         const std::vector<InOccurrence>& utab_rows,
                         const std::vector<InEdge>& etab_rows,
                         int k);
```

Logic:

1. Decide how many blocks this unitig needs based on Step 0's sizing.
2. Pack header bytes for each block.
3. 2-bit pack the DNA across blocks; the last 30 bytes (k-1 bp) of one block
   *do not need* to overlap into the next block — instead, the packer writes
   the unitig's bases sequentially across blocks, and the query path that
   extracts a k-mer near a block boundary reads across the boundary using the
   `next_block_offset` field.
4. Write UTAB rows into the head block.
5. Write ETAB rows into the head block. Resolve `to_block_offset` lazily — at
   pack time we don't know neighbor offsets yet (Step 5 fixes this).

**Deliverable**: a function that turns the Step 1 fixture into a printable
hex dump of bytes. Verify by hand on the smallest unitig.

---

## Step 5 — Flat Array assembler

**Goal**: glue all unitigs' blocks together into one global byte vector, and
back-patch the ETAB `to_block_offset` fields.

**Action**: Add `flat_assembler.hpp` / `.cpp`:

```cpp
struct FlatIndex {
    std::vector<std::byte>                    flat;          // the giant byte array
    std::unordered_map<uint32_t, uint32_t>    unitig_id_to_byte_offset; // build-only
    int                                        k;
};

FlatIndex assemble(const std::vector<InUnitig>& unitigs,
                   const std::vector<InOccurrence>& occs,
                   const std::vector<InEdge>& edges,
                   int k);
```

Logic:

1. Two-pass: first pass calls `pack_unitig` on every unitig with placeholder
   ETAB neighbor offsets, accumulates total byte size, records each unitig's
   start offset.
2. Second pass walks each unitig's ETAB rows and rewrites
   `to_block_offset = unitig_id_to_byte_offset[edge.to_id] / 64`.
3. After back-patching, drop the `unitig_id_to_byte_offset` map — the flat
   array is now self-describing.

**Deliverable**: in-memory FlatIndex object plus a `dump_flat_array(idx, path)`
helper that writes a hex dump for visual inspection.

---

## Step 6 — POS table + MPHF integration

**Goal**: build the canonical-k-mer → flat-array-offset map.

**Action**:

1. Iterate every unitig's bases, slide k-mer window, canonicalize each k-mer.
   For each canonical k-mer compute:
   - `target_block_byte_offset` (block within the flat array that contains
     this k-mer's start byte)
   - `local_offset_kmers` (k-mer index within that block's DNA)
   - assert `local_offset_kmers < 256` (the 8-bit constraint)

2. Pack as 64-bit POS entry:
   `pos_entry = (block_byte_offset << 8) | local_offset_kmers`

3. Build the MPHF over the canonical k-mers (use partner's `StaticMPHF` for
   now; replace with [BBHash](https://github.com/rizkg/BBHash) or
   [PTHash](https://github.com/jermp/pthash) later for production-scale
   keysets). Store `POS[mphf_slot] = pos_entry`.

**Deliverable**: a function `build_pos(flat_index, &mphf, &pos_table)` that
fills both. Verify by walking every k-mer of every unitig and checking
`POS[MPHF(canonical(K))]` recovers the right block + local offset.

---

## Step 7 — Query path

**Goal**: end-to-end "k-mer → list of (ref_id, ref_pos, orient)" lookup using
only the flat array + POS table + MPHF.

**Action**: Implement `std::vector<RefHit> query(K)`:

```text
1.  Kc = canonical(K)
2.  slot = MPHF(Kc)             ──► one cache line (MPHF data structure)
3.  pos  = POS[slot]              ──► one cache line (POS array)
4.  block_off = pos >> 8
5.  loc_kmers = pos & 0xFF
6.  read 64 bytes at flat[block_off]   ──► one cache line (the block!)
7.  verify: read_kmer(block.dna, loc_kmers·2 bits, k) == Kc
            (or RC(Kc) if orient flips)
8.  for each UTAB row in this block (in the same cache line):
      compute ref_pos via the formula:
        orient '+': ref_pos + (loc_kmers - entry_off)
        orient '-': ref_pos + (entry_off - loc_kmers)
      add to result list
9.  return result list
```

Step 8 is the only step where you need additional I/O if the unitig spans
multiple blocks and the UTAB rows live in the *head* block — in which case
follow `next_block_offset` *backwards* to the head. (For the prototype, store
`head_block_offset` in every block's header so step 8 needs at most one extra
fetch instead of a chain.)

**Deliverable**: `query(K)` returns the same answers as a brute-force scan over
the FASTA. Validate with a fixture-level test.

---

## Step 8 — Read-extension via ETAB

**Goal**: take a read longer than `k`, walk it through the index without
re-querying for every shifted k-mer.

**Action**: Implement `align_read(read)`:

```text
slot = MPHF(canonical(read[0..k]))    // initial seed
fetch block, position cursor at loc_kmers
for each next base b of the read:
    if cursor + 1 still inside this block's DNA:
        verify read[cursor..cursor+k] matches block DNA at cursor+1
        cursor++
    else:
        // unitig boundary
        for each ETAB row in block:
            if row.ext_base == b:
                jump to row.to_block_offset, reset cursor
                break
        else: report alignment break
```

**Deliverable**: `align_read("ACGT…")` returns the chain of unitig hops with
no hash-table lookups after the first. Within-unitig steps are pointer math.
Cross-unitig steps are exactly one cache-line fetch each.

---

## Step 9 — Validation harness

**Goal**: convince ourselves (and the reader of the report) the index is
correct before measuring performance.

**Action**: Add `validate.cpp`:

1. **Self-consistency**: for every k-mer of every unitig, query the index and
   assert the result includes the expected (block, local_off, orient) pair.
2. **Reference reconstruction**: walk each FASTA record k-mer by k-mer; the
   chain of UTAB hits should reproduce the original FASTA byte-for-byte.
3. **Fuzz**: generate 10k random k-mers, query both the flat index and a
   reference `std::unordered_map<string, vector<RefHit>>` brute-force index.
   Assert the result sets match.

**Deliverable**: `bin/validate` exits with `code 0` on a passing fixture,
prints the offending k-mer / position on failure.

---

## Step 10 — Benchmark harness

**Goal**: measure what the project promised — cache-miss reduction and query
throughput.

**Action**: Add `bench.cpp`:

1. **Build time**: `auto t0 = chrono::steady_clock::now(); build_index(...);`
   report `(t1 - t0)` for both the flat index and a baseline (e.g. partner's
   pre-flat MPHF + occurrence-list version, or a plain
   `unordered_map<string, vector<...>>`).
2. **Query throughput**: warm up the cache once, then measure
   `bench_query(idx, queries, n_iters)`. Report ns/query and queries/sec.
3. **Cache misses**: run under
   - Linux: `perf stat -e L1-dcache-load-misses,LLC-load-misses ./bench`
   - macOS: Instruments Counters template (or `xctrace record --template Counters`).
4. **Sweep**: parametrize by query count (10k, 100k, 1M) and by reference size,
   to show the design's scaling.

**Deliverable**: a `bench-report.md` with a table:

| index | build (s) | query (M ops/s) | L1 miss/q | LLC miss/q |
|---|---|---|---|---|
| baseline (partner's main.cpp) | … | … | … | … |
| flat index (this project)     | … | … | … | … |

These numbers go straight into the project report's "Results" section.

---

## Step 11 — Wire to the real pipeline (deferred until Steps 1–10 work)

Once everything works on the hardcoded fixture, swap the input source to TSVs
produced by `Project/bin/build_unitigs` and `Project/bin/build_occ`:

```bash
./bin/build_unitigs  200bp.fasta  data/200bp 31
./bin/build_occ      200bp.fasta  data/200bp 31
./bin/flatindex_build data/200bp                        # consumes the 4 TSVs
./bin/flatindex_query data/200bp.flat queries.fastq    # MPHF + POS + flat
```

Add TSV loaders that read the four files (`unitigs.tsv`, `edges.tsv`,
`occ.tsv`, `refs.tsv`) into the same `InUnitig` / `InOccurrence` / `InEdge`
structs the rest of the pipeline already speaks. Step 1's hardcoded fixture
remains useful as a regression test fixture forever.

---
STEP 12: Use Cuttlefish (https://github.com/COMBINE-lab/cuttlefish) to generate debrujin graph from a real genome and condurt rael performance testing...

