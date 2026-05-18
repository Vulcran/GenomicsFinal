#!/usr/bin/env python3
"""
Extract random read-length substrings from a FASTA file for benchmarking.

Produces a FASTQ file that works as input for both:
  - pufferfish align  -i <index> -1 queries.fastq
  - flatindex_query   <index.flat> --fastq queries.fastq

Reads are sampled from the actual reference sequence so every k-mer in
each read is guaranteed to be present in the index (realistic hit rate).
"""
import argparse
import random
import sys


def read_fasta(path):
    seqs = []
    current = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                if current:
                    seqs.append(''.join(current))
                    current = []
            elif line:
                current.append(line.upper())
    if current:
        seqs.append(''.join(current))
    return seqs


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('fasta', help='Input FASTA (e.g. ecoli.fasta)')
    p.add_argument('-n', '--num-queries', type=int, default=100_000,
                   help='Number of reads to generate (default: 100000)')
    p.add_argument('-l', '--read-len', type=int, default=100,
                   help='Read length in bp (default: 100; must be >= k)')
    p.add_argument('-o', '--output', default='queries.fastq',
                   help='Output FASTQ path (default: queries.fastq)')
    p.add_argument('--seed', type=int, default=42,
                   help='Random seed for reproducibility (default: 42)')
    args = p.parse_args()

    seqs = read_fasta(args.fasta)
    if not seqs:
        sys.exit('ERROR: no sequences found in FASTA')

    genome = ''.join(seqs)
    genome_len = len(genome)
    rl = args.read_len

    if genome_len < rl:
        sys.exit(f'ERROR: genome length {genome_len} < read length {rl}')

    valid = set('ACGT')
    rng = random.Random(args.seed)

    # Pre-scan: collect all positions where a window of length rl is clean.
    # For large genomes (E. coli ~4.6 MB) this is fast enough.
    print(f'Scanning {genome_len:,} bp for valid {rl}-bp windows...', file=sys.stderr)
    valid_starts = [
        i for i in range(genome_len - rl + 1)
        if all(b in valid for b in genome[i:i + rl])
    ]

    if not valid_starts:
        sys.exit(f'ERROR: no valid windows of length {rl} found (too many Ns?)')

    n = min(args.num_queries, len(valid_starts))
    if n < args.num_queries:
        print(f'Warning: only {len(valid_starts)} valid positions; '
              f'generating {n} queries', file=sys.stderr)

    chosen = rng.sample(valid_starts, n)
    qual = 'I' * rl  # Phred 40, arbitrary placeholder

    with open(args.output, 'w') as out:
        for i, pos in enumerate(chosen):
            seq = genome[pos:pos + rl]
            out.write(f'@read_{i} ref_pos={pos}\n{seq}\n+\n{qual}\n')

    print(f'Wrote {n:,} reads ({rl} bp) to {args.output}', file=sys.stderr)


if __name__ == '__main__':
    main()
