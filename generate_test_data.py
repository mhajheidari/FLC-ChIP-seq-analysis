#!/usr/bin/env python3
"""
generate_test_data.py

Generates a small, fully synthetic ChIP-seq dataset for CI smoke-testing:
- a tiny 2-chromosome genome FASTA
- a matching minimal GTF (a few genes, so ChIPseeker has something to
  annotate peaks against)
- paired-end FASTQ reads for every sample in the pipeline's sample sheet
  (P1-P4, Pi1-Pi4, C1, C2), split across two "sequencing run" directories
  to also exercise the multi-directory merge step

Three "peak" windows on chr1 have a strong artificial motif consensus
embedded in the genome sequence. Reads for the ChIP samples (P1-P4) are
enriched around these windows; input (Pi1-Pi4) and mock IP (C1, C2) reads
are uniform across the genome. This gives MACS3 genuine peaks to call and
STREME a genuine motif to discover, rather than just testing that the
scripts don't crash on empty/degenerate data.

This script has no dependencies beyond the Python standard library.
"""

import gzip
import os
import random
import argparse

RNG_SEED = 42

# Sample sheet: matches config.sh.example's names/ids
SAMPLES = {
    "P1": 346187, "P2": 346188, "P3": 346189, "P4": 346190,
    "Pi1": 346191, "Pi2": 346192, "Pi3": 346193, "Pi4": 346194,
    "C1": 346201, "C2": 346202,
}
CHIP_SAMPLES = ["P1", "P2", "P3", "P4"]

CHROM_LENGTHS = {"1": 8000, "2": 4000}
MOTIF = "AAAAATAGAAAGAAAAAAA"  # 19 bp artificial consensus, embedded at each peak center
PEAK_WINDOW = 150
PEAK_CENTERS = [1000, 3000, 5500]  # all on chr1, well clear of chromosome ends

READ_LENGTH = 75
FRAGMENT_LENGTH = 200
N_FRAGMENTS_PER_SAMPLE = 3000
ENRICHED_FRACTION = 0.85  # fraction of ChIP-sample fragments drawn from peak windows
PEAK_SPREAD = 40  # stdev (bp) of fragment center around each peak center, for ChIP samples


COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def random_seq(length, rng):
    return "".join(rng.choice("ACGT") for _ in range(length))


def build_genome(rng):
    genome = {}
    for chrom, length in CHROM_LENGTHS.items():
        genome[chrom] = list(random_seq(length, rng))

    # Embed the motif at each peak center on chr1
    for center in PEAK_CENTERS:
        start = center - len(MOTIF) // 2
        genome["1"][start:start + len(MOTIF)] = list(MOTIF)

    return {chrom: "".join(seq) for chrom, seq in genome.items()}


def write_fasta(genome, out_path, line_width=60):
    with open(out_path, "w") as f:
        for chrom, seq in genome.items():
            f.write(f">{chrom}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i + line_width] + "\n")


def write_gtf(out_path):
    # A few minimal genes on chr1 (overlapping/near the peak windows, so
    # ChIPseeker has real promoter/genic context to annotate against) and
    # one on chr2 (so both chromosomes have at least one feature).
    genes = [
        ("1", 800, 1400, "+", "gene1"),
        ("1", 2700, 3400, "-", "gene2"),
        ("1", 5300, 6000, "+", "gene3"),
        ("2", 1000, 2000, "+", "gene4"),
    ]
    with open(out_path, "w") as f:
        for chrom, start, end, strand, gene_id in genes:
            attrs_gene = f'gene_id "{gene_id}"; gene_name "{gene_id}";'
            attrs_tx = f'gene_id "{gene_id}"; transcript_id "{gene_id}.1";'
            f.write(f"{chrom}\ttest\tgene\t{start}\t{end}\t.\t{strand}\t.\t{attrs_gene}\n")
            f.write(f"{chrom}\ttest\ttranscript\t{start}\t{end}\t.\t{strand}\t.\t{attrs_tx}\n")
            f.write(f"{chrom}\ttest\texon\t{start}\t{end}\t.\t{strand}\t.\t{attrs_tx}\n")


def sample_fragment_start(sample_name, genome_lengths, rng):
    if sample_name in CHIP_SAMPLES and rng.random() < ENRICHED_FRACTION:
        center = rng.choice(PEAK_CENTERS)
        pos = int(rng.gauss(center - FRAGMENT_LENGTH // 2, PEAK_SPREAD))
        pos = max(10, min(pos, genome_lengths["1"] - FRAGMENT_LENGTH - 10))
        return "1", pos
    else:
        # Uniform background: pick a chromosome weighted by how much valid
        # fragment-start range it has, then a uniform position within it.
        valid_ranges = {c: max(0, L - FRAGMENT_LENGTH) for c, L in genome_lengths.items()}
        total = sum(valid_ranges.values())
        r = rng.uniform(0, total)
        cum = 0
        for c, v in valid_ranges.items():
            cum += v
            if r <= cum:
                return c, (rng.randrange(0, v) if v > 0 else 0)
        # Should be unreachable, but fail safe rather than crash
        return "1", 10


def mutate(seq, rate, rng):
    bases = "ACGT"
    seq = list(seq)
    for i in range(len(seq)):
        if rng.random() < rate:
            seq[i] = rng.choice(bases)
    return "".join(seq)


def simulate_reads(sample_name, genome, rng, n_fragments=N_FRAGMENTS_PER_SAMPLE):
    r1_records = []
    r2_records = []
    genome_lengths = {c: len(s) for c, s in genome.items()}

    for i in range(n_fragments):
        chrom, start = sample_fragment_start(sample_name, genome_lengths, rng)
        frag = genome[chrom][start:start + FRAGMENT_LENGTH]

        r1_seq = mutate(frag[:READ_LENGTH], 0.01, rng)
        r2_seq = mutate(revcomp(frag[-READ_LENGTH:]), 0.01, rng)

        qual = "I" * READ_LENGTH
        read_name = f"{sample_name}:{chrom}:{start}:{i}"
        r1_records.append(f"@{read_name} 1\n{r1_seq}\n+\n{qual}\n")
        r2_records.append(f"@{read_name} 2\n{r2_seq}\n+\n{qual}\n")

    return r1_records, r2_records


def write_fastq_gz(records, out_path):
    with gzip.open(out_path, "wt") as f:
        f.writelines(records)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="test_data", help="Output directory")
    args = ap.parse_args()

    rng = random.Random(RNG_SEED)

    run1_dir = os.path.join(args.outdir, "sequencing_run_1")
    run2_dir = os.path.join(args.outdir, "sequencing_run_2")
    os.makedirs(run1_dir, exist_ok=True)
    os.makedirs(run2_dir, exist_ok=True)

    genome = build_genome(rng)
    write_fasta(genome, os.path.join(args.outdir, "test_genome.fa"))
    write_gtf(os.path.join(args.outdir, "test_genome.gtf"))

    total_genome_size = sum(len(s) for s in genome.values())

    for sample_name, sample_id in SAMPLES.items():
        r1, r2 = simulate_reads(sample_name, genome, rng)

        half = len(r1) // 2
        write_fastq_gz(r1[:half], os.path.join(
            run1_dir, f"{sample_name}_S1_{sample_id}_L001_R1_001.fastq.gz"))
        write_fastq_gz(r2[:half], os.path.join(
            run1_dir, f"{sample_name}_S1_{sample_id}_L001_R2_001.fastq.gz"))
        write_fastq_gz(r1[half:], os.path.join(
            run2_dir, f"{sample_name}_S1_{sample_id}_L002_R1_001.fastq.gz"))
        write_fastq_gz(r2[half:], os.path.join(
            run2_dir, f"{sample_name}_S1_{sample_id}_L002_R2_001.fastq.gz"))

        print(f"  {sample_name} (id {sample_id}): {len(r1)} read pairs written")

    print(f"\nGenome: {args.outdir}/test_genome.fa ({total_genome_size} bp total)")
    print(f"GTF:    {args.outdir}/test_genome.gtf")
    print(f"Effective genome size for MACS3 -g: {total_genome_size}")


if __name__ == "__main__":
    main()
