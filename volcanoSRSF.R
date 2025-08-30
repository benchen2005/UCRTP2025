###############################################################################
# Volcano – highlight SRSF family (robust matching, large fonts)
###############################################################################

# --- User paths / thresholds ---
file_path  <- "/Users/ben/Documents/UCRTP/Filtered Events/noFilterFDR.SE.tsv"
psi_thresh <- 0.10
fdr_thresh <- 0.05
out_png    <- "volcano_SRSF_family.png"
out_pdf    <- "volcano_SRSF_family.pdf"

# --- Packages (namespaced to avoid conflicts) ---
pkgs <- c("ggplot2","ggrepel","readr","dplyr")
new  <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new)
invisible(lapply(pkgs, library, character.only = TRUE))

numify <- function(x) if (is.character(x)) readr::parse_number(x) else as.numeric(x)

# --- Load & clean ---
stopifnot(file.exists(file_path))
se <- utils::read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)

stopifnot(all(c("geneSymbol","PSI","FDR") %in% names(se)))
se <- se |>
  dplyr::mutate(
    geneSymbol = trimws(geneSymbol),
    PSI        = numify(PSI),
    FDR        = numify(FDR)
  )

# handle FDR==0 then keep only finite rows
if (any(se$FDR == 0, na.rm = TRUE)) {
  eps <- min(se$FDR[se$FDR > 0], na.rm = TRUE) * 0.5
  se$FDR[se$FDR == 0] <- eps
}
se <- se |>
  dplyr::filter(is.finite(PSI), is.finite(FDR), FDR >= 0)

# --- Robust SRSF detection (strict + loose + aliases) ---
syms     <- toupper(se$geneSymbol)
strict   <- grepl("^SRSF[0-9]+$", syms)              # exact SRSF#
loose    <- grepl("SRSF", syms)                      # contains SRSF anywhere
aliases1 <- syms %in% c("SF2","ASF","ASF/SF2")       # → SRSF1
aliases2 <- syms %in% c("SC35","SRP30B")             # → SRSF2

# canonical label
label <- rep(NA_character_, length(syms))
label[strict]  <- syms[strict]
label[aliases1] <- "SRSF1"
label[aliases2] <- "SRSF2"

isSRSF <- strict | aliases1 | aliases2 | loose

# Prefer canonical label where known; otherwise keep original symbol for SRSF rows
label[is.na(label) & isSRSF] <- se$geneSymbol[is.na(label) & isSRSF]
se$isSRSF <- isSRSF
se$label  <- label

# quick sanity output
cat("Detected SRSF rows:", sum(se$isSRSF), " / ", nrow(se), "\n")
print(utils::head(sort(unique(se$label[se$isSRSF]))))

# --- Y zoom without dropping points ---
y_cap <- stats::quantile(-log10(se$FDR), probs = 0.995, na.rm = TRUE) * 1.02

# --- Plot (bigger fonts & labels) ---
p <- ggplot2::ggplot(se, ggplot2::aes(x = PSI, y = -log10(FDR))) +
  # background points
  ggplot2::geom_point(data = subset(se, !isSRSF),
                      colour = "grey80", size = 1.4, alpha = 0.7) +
  # SRSF highlights
  ggplot2::geom_point(data = subset(se, isSRSF),
                      colour = "dodgerblue3", fill = "dodgerblue2",
                      shape = 21, size = 3.6, stroke = 0.5) +
  ggrepel::geom_text_repel(
    data = subset(se, isSRSF),
    ggplot2::aes(label = label),
    size = 6.5,                      # larger text
    fontface = "bold",
    min.segment.length = 0,
    box.padding  = 0.3,
    point.padding = 0.25,
    max.overlaps = Inf
  ) +
  # thresholds
  ggplot2::geom_vline(xintercept = c(-psi_thresh, psi_thresh),
                      linetype = "dashed", linewidth = 0.6) +
  ggplot2::geom_hline(yintercept = -log10(fdr_thresh),
                      linetype = "dashed", linewidth = 0.6) +
  # zoom
  ggplot2::coord_cartesian(ylim = c(0, y_cap)) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.07))) +
  # labels/theme (bigger base font)
  ggplot2::labs(
    title    = "Volcano plot – SRSF family splicing events",
    subtitle = "SRSF1–SRSF7 & SRSF10 highlighted (aliases: SF2/ASF→SRSF1; SC35/SRp30b→SRSF2)",
    x        = "PSI difference (ΔPSI)",
    y        = expression(-log[10]~FDR)
  ) +
  ggplot2::theme_minimal(base_size = 18) +
  ggplot2::theme(
    plot.title    = ggplot2::element_text(face = "bold", size = 20),
    plot.subtitle = ggplot2::element_text(size = 16),
    panel.grid.minor = ggplot2::element_blank()
  )

print(p)
ggplot2::ggsave(out_png, p, width = 12, height = 9, dpi = 300)
ggplot2::ggsave(out_pdf, p, width = 12, height = 9)
cat("Saved:", normalizePath(out_png), "\n")
