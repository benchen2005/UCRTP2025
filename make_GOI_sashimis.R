###############################################################################
# 0 ── PATHS                                                                  #
###############################################################################
root_dir  <- "/Users/ben/Documents/UCRTP/Dutta Drive"
bam_dir   <- file.path(root_dir, "bam")
tsv_file  <- file.path(root_dir, "tsv/Filtered SE copy.tsv")
gtf_file  <- file.path(root_dir, "gtf/gencode.v48.primary_assembly.annotation.gtf")
out_dir   <- file.path(root_dir, "plots/sashimi_GOI")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 1 ── GENES OF INTEREST                                                      #
###############################################################################
goi <- c("SRSF5","SRSF2","SRSF6","SRSF3","SRSF11","BRCA1","RSRP1","Sirt2","Sirt6",
         "NAMPT","HNRNPA1","BCL2L1","BCL2","SLC7A11","UBE2F","ESRP1","HNRNPD",
         "HNRNPU","UBE3A","UBE4B","NOTCH1","SIRT7","HRAS","PFKL","XBP1","NAFKB2",
         "TAB3","TGFB","MAT2A","RSP11","RTN2","STK36","MLXIPL","CPT2","KMT2C",
         "MALT1","HNRNPK","EZH2","PTPN11","CTNNB1","CRTC2","EZH1","KAT2A",
         "NDUFS2","YTHDF2","LDAH","HDAC6","NDUFAF6","BTAF1","RBM39","ACSS2",
         "UPF3A","NCOR1","ATM","YAP1","MALAT1","ME2")

###############################################################################
# 2 ── LOAD PACKAGES                                                          #
###############################################################################
pkgs <- c("dplyr","maser","Rsamtools")
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(need, ask = FALSE)
}
lapply(pkgs, library, character.only = TRUE)

###############################################################################
# 3 ── READ SE TSV                                                            #
###############################################################################
stopifnot(file.exists(tsv_file))
se <- read.delim(tsv_file, check.names = FALSE, stringsAsFactors = FALSE)

# build ID if missing
if (!"ID" %in% names(se)) {
  c_chr   <- grep("^chr$|chrom", names(se), value = TRUE)[1]
  c_start <- grep("start", names(se), value = TRUE)[1]
  c_end   <- grep("end",   names(se), value = TRUE)[1]
  se$ID   <- with(se, paste0(geneSymbol, "_", .data[[c_chr]], ":",
                             .data[[c_start]], "-", .data[[c_end]]))
}

###############################################################################
# 4 ── EFFECT-SIZE COLUMN (PSI)                                               #
###############################################################################
abs_col <- "PSI"
stopifnot(abs_col %in% names(se))

###############################################################################
# 5 ── SELECT ONE EVENT PER GOI                                               #
###############################################################################
goi_events <- se %>%
  filter(geneSymbol %in% goi) %>%
  mutate(absPSI = abs(as.numeric(.data[[abs_col]]))) %>%
  group_by(geneSymbol) %>%
  slice_max(absPSI, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("▶ Will generate", nrow(goi_events), "sashimi plots (1 per gene)\n")

###############################################################################
# 6 ── PREP BAM FILES & INDEX IF NEEDED                                       #
###############################################################################
bam_files <- list.files(bam_dir, "\\.bam$", full.names = TRUE)
stopifnot(length(bam_files) > 0)

for (bf in bam_files)
  if (!file.exists(paste0(bf, ".bai"))) Rsamtools::indexBam(bf)

sample_labels <- factor(sub("\\.bam$","",basename(bam_files)))
o <- order(sample_labels)
bam_files     <- bam_files[o]
sample_labels <- sample_labels[o]

###############################################################################
# 7 ── LOOP & PLOT SASHIMIS                                                   #
###############################################################################
for (i in seq_len(nrow(goi_events))) {
  ev   <- goi_events$ID[i]
  gene <- goi_events$geneSymbol[i]
  out_pdf <- file.path(out_dir, paste0(gene, "_", ev, ".pdf"))
  cat(sprintf("  [%d/%d] %s\n", i, nrow(goi_events), gene))
  tryCatch(
    maser::plotSashimi(
      event      = ev,
      data       = se,
      bamfiles   = bam_files,
      labels     = sample_labels,
      gtfFile    = gtf_file,
      minDepth   = 5,
      outFile    = out_pdf
    ),
    error = function(e) message("    skipped: ", e$message)
  )
}

cat("\n✔ PDFs saved in", normalizePath(out_dir), "\n")
