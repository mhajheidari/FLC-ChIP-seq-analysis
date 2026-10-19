# ============================================================
# Find reproducible ChIP-seq peaks for the FLC::FLCg-mCherry ChIP
# samples (P1-P4, pooled biological replicates), apply the mock-IP
# specificity filter, and annotate both peak sets with ChIPseeker
# (Arabidopsis thaliana TAIR10)
#
# Reproducibility definition (see manuscript Methods):
#   A peak is retained only if it is called (MACS3, q < 0.05) in
#   ALL FOUR pooled biological replicates (P1, P2, P3, P4). Peaks from
#   the four replicates are merged into a union of non-overlapping
#   regions (GenomicRanges::reduce), and a region is kept only if it
#   overlaps a called peak from every one of the four replicates.
#
# Specificity filter:
#   Regions that also overlap a peak called in the mock IP samples
#   (C1/C2, RFP-Trap on non-FLCg-mCherry chromatin) are removed as
#   likely nonspecific RFP-Trap/bead-associated signal.
#
# Requires run_ChIPseq_pipeline.sh to have already been run (it calls
# peaks for both the P1-P4 vs Pi1-Pi4 pairs and the mock IP samples
# C1/C2 vs the pooled input background).
# ============================================================

library(ChIPseeker)
library(GenomicRanges)
library(GenomeInfoDb)
library(txdbmaker)
library(AnnotationDbi)
library(rtracklayer)

source("config.R")

PEAKDIR <- file.path(OUTDIR, "reproducible_peaks")
ANNOTDIR <- file.path(OUTDIR, "peak_annotation_reproducible")

dir.create(PEAKDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ANNOTDIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# PART 1: FIND REPRODUCIBLE, MOCK-FILTERED PEAKS
# ============================================================

# ------------------------------------------------------------
# Helper: read a MACS3 narrowPeak file into a GRanges
# ------------------------------------------------------------
read_narrowpeak <- function(sample_name) {
    peak_file <- file.path(
        OUTDIR, "macs3", sample_name,
        paste0(sample_name, "_peaks.narrowPeak")
    )

    if (!file.exists(peak_file)) {
        stop("Missing MACS3 output for ", sample_name, ": ", peak_file)
    }

    # MACS3 writes an empty narrowPeak file when no peaks pass the cutoff
    # (e.g. a clean mock IP); read.table() would error on that.
    if (file.size(peak_file) == 0) {
        return(GRanges(
            seqnames = character(0),
            ranges = IRanges(start = integer(0), end = integer(0)),
            qvalue = numeric(0)
        ))
    }

    df <- read.table(peak_file, header = FALSE, stringsAsFactors = FALSE)
    colnames(df) <- c(
        "chrom", "start", "end", "name", "score", "strand",
        "signalValue", "pValue", "qValue", "summit"
    )

    GRanges(
        seqnames = df$chrom,
        ranges = IRanges(start = df$start + 1, end = df$end),  # BED 0-based -> 1-based
        qvalue = df$qValue
    )
}

write_bed4 <- function(gr, out_file) {
    df <- data.frame(
        chrom = as.character(seqnames(gr)),
        start = start(gr) - 1L,  # back to BED 0-based
        end = end(gr),
        qvalue = gr$qvalue
    )
    write.table(
        df, file = out_file, sep = "\t",
        quote = FALSE, row.names = FALSE, col.names = FALSE
    )
}

cat("=== Reading MACS3 peaks for", paste(CHIP_REPLICATES, collapse = ", "), "===\n")

rep_peaks <- lapply(CHIP_REPLICATES, read_narrowpeak)
names(rep_peaks) <- CHIP_REPLICATES

for (r in CHIP_REPLICATES) {
    cat("  ", r, ":", length(rep_peaks[[r]]), "peaks\n")
}

cat("=== Computing reproducible peaks (present in all",
    length(CHIP_REPLICATES), "replicates) ===\n")

all_peaks <- do.call(c, unname(rep_peaks))
union_regions <- reduce(all_peaks)

keep <- rep(TRUE, length(union_regions))
for (r in CHIP_REPLICATES) {
    keep <- keep & overlapsAny(union_regions, rep_peaks[[r]])
}

reproducible_peaks <- union_regions[keep]

# Representative q-value per reproducible region: the least
# significant (maximum) q-value among the overlapping per-replicate
# peaks, i.e. a conservative combined estimate.
rep_qvalue <- sapply(seq_along(reproducible_peaks), function(i) {
    region <- reproducible_peaks[i]
    qvals <- sapply(CHIP_REPLICATES, function(r) {
        hits <- findOverlaps(region, rep_peaks[[r]])
        max(rep_peaks[[r]]$qvalue[subjectHits(hits)])
    })
    max(qvals)
})
reproducible_peaks$qvalue <- rep_qvalue

cat("Reproducible peaks (before mock filtering):", length(reproducible_peaks), "\n")

write_bed4(reproducible_peaks, file.path(PEAKDIR, "P_reproducible_peaks.bed"))
cat("Written:", file.path(PEAKDIR, "P_reproducible_peaks.bed"), "\n\n")

cat("=== Reading mock IP peaks (", paste(MOCK_SAMPLES, collapse = ", "), ") ===\n")

mock_peaks_list <- lapply(MOCK_SAMPLES, read_narrowpeak)
mock_peaks <- reduce(do.call(c, mock_peaks_list))
cat("Mock IP peaks (union of", paste(MOCK_SAMPLES, collapse = "+"), "):",
    length(mock_peaks), "\n")

nonspecific <- overlapsAny(reproducible_peaks, mock_peaks)
mock_filtered_peaks <- reproducible_peaks[!nonspecific]

cat("Reproducible peaks removed as nonspecific (overlap mock IP):",
    sum(nonspecific), "\n")
cat("Reproducible peaks after mock filtering:",
    length(mock_filtered_peaks), "\n")

write_bed4(
    mock_filtered_peaks,
    file.path(PEAKDIR, "P_reproducible_peaks_mockfiltered.bed")
)
cat("Written:", file.path(PEAKDIR, "P_reproducible_peaks_mockfiltered.bed"), "\n\n")

# ============================================================
# PART 2: ANNOTATE BOTH PEAK SETS WITH CHIPSEEKER
# ============================================================

cat("=== Creating TAIR10 TxDb ===\n")

txdb <- txdbmaker::makeTxDbFromGFF(GTF, format = "gtf")

cat("TxDb created successfully.\n")
cat("Chromosomes in TxDb:\n")
print(seqlevels(txdb))

txdb_file <- file.path(ANNOTDIR, "TAIR10.sqlite")
AnnotationDbi::saveDb(txdb, file = txdb_file)
cat("TxDb saved to:", txdb_file, "\n\n")

# Draw a plot to a PDF; on failure, warn and remove the partial file
# instead of stopping the script.
save_plot <- function(file, plot_fun, width, height) {
    pdf(file, width = width, height = height)
    ok <- tryCatch({
        print(plot_fun())
        TRUE
    }, error = function(e) {
        warning("Could not draw ", basename(file), ": ", conditionMessage(e),
                call. = FALSE)
        FALSE
    })
    dev.off()
    if (!ok && file.exists(file)) file.remove(file)
    invisible(ok)
}

peak_sets <- c(
    "P_reproducible_peaks",
    "P_reproducible_peaks_mockfiltered"
)

for (set_name in peak_sets) {

    cat("\n============================================\n")
    cat("Processing:", set_name, "\n")
    cat("============================================\n")

    peak_file <- file.path(PEAKDIR, paste0(set_name, ".bed"))

    if (!file.exists(peak_file)) {
        warning("Peak file does not exist: ", peak_file)
        next
    }

    df <- read.table(peak_file, header = FALSE, stringsAsFactors = FALSE)
    colnames(df) <- c("chr", "start", "end", "qvalue")

    peaks <- GRanges(
        seqnames = df$chr,
        ranges = IRanges(start = df$start + 1, end = df$end),   # BED 0-based to 1-based
        qvalue = df$qvalue
    )

    cat("Number of peaks:", length(peaks), "\n")

    common_chr <- intersect(seqlevels(peaks), seqlevels(txdb))
    cat("Chromosomes shared with TxDb:", paste(common_chr, collapse = ", "), "\n")

    # Promoter window (TSS_REGION) and other annotation parameters come
    # from config.R
    peakAnno <- annotatePeak(
        peaks,
        TxDb = txdb,
        tssRegion = TSS_REGION,
        level = "gene",
        assignGenomicAnnotation = TRUE,
        genomicAnnotationPriority = c(
            "Promoter", "5UTR", "3UTR", "Exon", "Intron",
            "Downstream", "Intergenic"
        ),
        addFlankGeneInfo = TRUE,
        flankDistance = 5000,
        sameStrand = FALSE,
        verbose = TRUE
    )

    anno_df <- as.data.frame(peakAnno)

    outfile <- file.path(ANNOTDIR, paste0(set_name, "_annotation.tsv"))
    write.table(anno_df, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Annotation table written:\n")
    cat(outfile, "\n")

    # The plots are a convenience; the annotation table above is the real
    # result. A plot that cannot be drawn (e.g. ChIPseeker's distance-to-TSS
    # plot needs peaks on both sides of a TSS) must not abort the pipeline.
    save_plot(
        file.path(ANNOTDIR, paste0(set_name, "_annotation_distribution.pdf")),
        function() plotAnnoBar(peakAnno), width = 8, height = 6
    )
    save_plot(
        file.path(ANNOTDIR, paste0(set_name, "_distance_to_TSS.pdf")),
        function() plotDistToTSS(peakAnno), width = 8, height = 6
    )
    save_plot(
        file.path(ANNOTDIR, paste0(set_name, "_annotation_pie.pdf")),
        function() plotAnnoPie(peakAnno), width = 7, height = 7
    )
}

cat("\n============================================\n")
cat("REPRODUCIBLE PEAK FINDING + ANNOTATION FINISHED\n")
cat("============================================\n")
cat("Final peak set for motif analysis:",
    file.path(PEAKDIR, "P_reproducible_peaks_mockfiltered.bed"), "\n")
