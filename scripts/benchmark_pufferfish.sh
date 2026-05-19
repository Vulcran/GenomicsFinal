#!/usr/bin/env bash
# benchmark_pufferfish.sh — time flat index vs real Pufferfish on the same data.
#
# Run from the GenomicsFinal/ project root on macOS (Apple Silicon or Intel):
#
#   bash scripts/benchmark_pufferfish.sh ecoli.fasta [k] [n_queries]
#
# Positional arguments (all optional):
#   $1  Path to reference FASTA           (default: ecoli.fasta)
#   $2  k-mer size                        (default: 31)
#   $3  Number of query reads to generate (default: 100000)
#
# Outputs timing for:
#   - flatindex_query --fastq  (your flat index)
#   - pufferfish align         (original Pufferfish)
#
# Requirements: cmake, python3, conda (for pufferfish install)
# Pufferfish is installed automatically if not found.

set -euo pipefail

# ── Config ────────────────────────────────────────────────────────────────────
FASTA="${1:-ecoli.fasta}"
K="${2:-31}"
N_QUERIES="${3:-100000}"
READ_LEN=100
DATA_DIR="data"
PREFIX="${DATA_DIR}/ecoli"
QUERIES="${DATA_DIR}/queries.fastq"
KMERS="${DATA_DIR}/queries.kmers"
FLAT_INDEX="${DATA_DIR}/ecoli.flat"
PUF_INDEX="${DATA_DIR}/puf_index"
THREADS="$(sysctl -n hw.logicalcpu)"

die() { echo "ERROR: $*" >&2; exit 1; }
hr()  { echo "──────────────────────────────────────────────────────────"; }

hr
echo " Flat Index vs Pufferfish Benchmark"
printf "  Reference : %s\n" "$FASTA"
printf "  k         : %s\n" "$K"
printf "  Queries   : %s reads x %s bp\n" "$N_QUERIES" "$READ_LEN"
printf "  Build threads : %s (query comparison is single-threaded)\n" "$THREADS"
hr
echo

[ -f "$FASTA" ] || die "FASTA not found: '$FASTA'  (pass full path as first argument)"
mkdir -p "$DATA_DIR"

# ── 1. Build the project ──────────────────────────────────────────────────────
echo "[1/6] Building project (Release)..."
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -Wno-dev -DCMAKE_EXPORT_COMPILE_COMMANDS=OFF \
    > /dev/null 2>&1
cmake --build build --config Release -j"$THREADS" 2>&1 \
    | grep -E 'error:|Linking|Building' || true
[ -f "build/flatindex_query" ] || [ -f "build/Release/flatindex_query.exe" ] \
    || die "Build failed — flatindex_query not found"
# Normalise binary path (macOS vs Windows)
FIQ="build/flatindex_query"
FIB="build/flatindex_build"
BU="build/build_unitigs"
BO="build/build_occ"
echo "    ok"
echo

# ── 2. Generate queries ───────────────────────────────────────────────────────
echo "[2/6] Generating $N_QUERIES x ${READ_LEN}-bp queries from reference..."
python3 scripts/generate_queries.py \
    "$FASTA" \
    --num-queries "$N_QUERIES" \
    --read-len    "$READ_LEN" \
    --output      "$QUERIES"
echo

# ── 3. Build flat index ───────────────────────────────────────────────────────
echo "[3/6] Building flat index..."
"$BU" "$FASTA" "$PREFIX" "$K" 2>&1 | grep '\[build_unitigs\]' | tail -2
"$BO" "$FASTA" "$PREFIX" "$K" 2>&1 | grep -E 'sanity|reconstruction|wrote' || true

echo "── Flat index build time (flatindex_build only) ─────"
{ time "$FIB" --tsv "$PREFIX" --k "$K" --out "$FLAT_INDEX" 2>&1 | tail -2; } 2>&1
echo

# ── 4. Install pufferfish if needed ──────────────────────────────────────────
echo "[4/6] Checking pufferfish..."
if ! command -v pufferfish &>/dev/null; then
    echo "    not found — installing via conda..."

    # Try native arm64 first (works on newer bioconda builds)
    if conda install -y -c bioconda -c conda-forge pufferfish 2>/dev/null; then
        echo "    installed (arm64 native)"
    else
        echo "    arm64 unavailable; falling back to osx-64 (Rosetta 2)..."
        # Rosetta 2 must be installed: softwareupdate --install-rosetta
        CONDA_SUBDIR=osx-64 conda install -y -c bioconda -c conda-forge pufferfish \
            || die "conda install failed.
  Manual build option:
    brew install cmake boost jellyfish
    git clone https://github.com/COMBINE-lab/pufferfish
    cd pufferfish && mkdir build && cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make -j$THREADS
    export PATH=\$PWD:\$PATH"
        echo "    installed (osx-64 via Rosetta)"
    fi
else
    echo "    found: $(command -v pufferfish)"
fi
echo

# ── 5. Build pufferfish index ─────────────────────────────────────────────────
echo "[5/6] Building pufferfish index (k=$K)..."
echo "── Pufferfish index build time ──────────────────────"
mkdir -p "$PUF_INDEX"
{ time pufferfish index \
    -r "$FASTA" \
    -o "$PUF_INDEX" \
    -k "$K" \
    -p 1 \
    > /dev/null 2>/dev/null; } 2>&1
echo

# ── 6. Timed comparison ───────────────────────────────────────────────────────
echo "[6/6] Timed query comparison ($N_QUERIES reads)..."
echo

# Flat index — align_read() walk for each read in the FASTQ
echo "── Your flat index (flatindex_query --fastq) ──"
{ time "$FIQ" "$FLAT_INDEX" --fastq "$QUERIES" --noOutput; } 2>&1
echo

# Pufferfish — seed-and-extend alignment, SAM output suppressed
echo "── Original Pufferfish (pufferfish align) ──────"
{ time pufferfish align \
    -i "$PUF_INDEX" \
    --read "$QUERIES" \
    -t 1 \
    --noOutput \
    > /dev/null 2>/dev/null; } 2>&1
echo

hr
echo " Alignment comparison done."
hr
echo

# ── 7. K-mer lookup comparison ────────────────────────────────────────────────
echo "[7/7] K-mer lookup comparison ($N_QUERIES individual k-mers)..."
echo

# Extract first k-mer from each read in the FASTQ
awk 'NR%4==2 {print substr($0,1,'"$K"')}' "$QUERIES" > "$KMERS"

echo "── Your flat index (flatindex_query --kmer-file) ──"
{ time "$FIQ" "$FLAT_INDEX" --kmer-file "$KMERS" --noOutput; } 2>&1
echo

echo "── Original Pufferfish (pufferfish kquery) ─────────"
{ time pufferfish kquery \
    -i "$PUF_INDEX" \
    -q "$KMERS" \
    -p 1 \
    > /dev/null 2>/dev/null; } 2>&1
echo

hr
echo " Copy the 'real' wall-clock times above into bench-report.md."
echo " Both tools queried the same $N_QUERIES k-mers with k=$K."
hr
