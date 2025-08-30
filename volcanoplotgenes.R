###############################################################################
# 0 ── user settings                                                          #
###############################################################################
file_path   <- "/Users/ben/Documents/UCRTP/Filtered Events/Filtered SE.tsv"   # <- path to TSV
psi_thresh  <- 0.10   # |ΔPSI| threshold
fdr_thresh  <- 0.001  # *** updated ***
genes_interest <- c(
  "SRSF5","SRSF2","SRSF4","SRSF6","SRSF3","SRSF11","BRCA1","RSRP1","Sirt2","Sirt6",
  "NAMPT","HNRNPA1","BCL2L1","BCL2","SLC7A11","UBE2F","ESRP1","HNRNPD",
  "HNRNPU","UBE3A","UBE4B","NOTCH1","SIRT7","HRAS","PFKL","XBP1"
)

###############################################################################
# 1 ── load packages (install if missing)                                     #
###############################################################################
pkgs <- c("ggplot2","ggrepel","readr","dplyr")
new  <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new)) install.packages(new)
invisible(lapply(pkgs, library, character.only = TRUE))

###############################################################################
# 2 ── read & clean                                                           #
###############################################################################
stopifnot(file.exists(file_path))

se <- read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)
req_cols <- c("geneSymbol","PSI","FDR")
stopifnot(all(req_cols %in% names(se)))

clean_num <- function(x) if (is.character(x)) readr::parse_number(x) else as.numeric(x)
se <- se %>%
  mutate(
    PSI = clean_num(PSI),
    FDR = clean_num(FDR)
  )

###############################################################################
# 3 ── flag and subset                                                        #
###############################################################################
se <- se %>%
  mutate(
    Sig   = abs(PSI) > psi_thresh & FDR < fdr_thresh,
    GOI   = geneSymbol %in% genes_interest,              # Gene-of-interest?
    label = ifelse(GOI, geneSymbol, NA)
  )

###############################################################################
# 4 ── build the volcano plot                                                 #
###############################################################################
base_plot <- ggplot(se, aes(PSI, -log10(FDR))) +
  # non-significant layer first
  geom_point(data = subset(se, !Sig),
             colour = "grey75", size = 1) +
  # significant but NOT GOI
  geom_point(data = subset(se, Sig & !GOI),
             colour = "#F8766D", size = 1.2) +
  # highlight genes-of-interest
  geom_point(data = subset(se, GOI),
             colour = "dodgerblue3", fill = "dodgerblue2",
             size = 3, shape = 21, stroke = 0.4) +
  geom_text_repel(data = subset(se, GOI),
                  aes(label = label),
                  min.segment.length = 0,
                  size = 3,
                  box.padding  = 0.3,
                  point.padding = 0.15,
                  max.overlaps = Inf,
                  fontface = "bold") +
  # threshold lines
  geom_vline(xintercept = c(-psi_thresh, psi_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_thresh),           linetype = "dashed") +
  # axis / theme
  labs(title = "Volcano plot – Skipped-exon (SE) events",
       subtitle = paste("Red = Sig (|ΔPSI| >", psi_thresh,
                        "& FDR <", fdr_thresh, ")  •  Blue = Genes of interest"),
       x = "PSI difference (ΔPSI)",
       y = expression(-log[10]~FDR)) +
  theme_minimal(base_size = 11)

print(base_plot)

###############################################################################
# 5 ── save hi-res output                                                     #
###############################################################################
ggsave("volcano_SE_GOI.png", base_plot, width = 9, height = 9, dpi = 300)
ggsave("volcano_SE_GOI.pdf", base_plot, width = 9, height = 9)

###############################################################################
# 6 ── quick interpretation/guide (prints to console)                         #
###############################################################################
cat("\nINTERPRETATION\n",
    "• Each dot is a skipped-exon event.\n",
    "• X-axis (ΔPSI): negative = exon included LESS in knock-down than control;\n",
    "                 positive = exon included MORE.\n",
    "• Y-axis: higher = more statistically significant (smaller FDR).\n",
    "• Dashed vertical lines at ±", psi_thresh,
    " show the effect-size threshold.\n",
    "• Dashed horizontal line at –log10(", fdr_thresh,
    ") ≈ ", round(-log10(fdr_thresh), 2),
    " shows the FDR cut-off.\n",
    "• Grey = not significant; light-red = significant; blue = your genes of\n",
    "  interest (also significant if right of/below thresholds). Labels use\n",
    "  bold text and repel to avoid overlap.\n", sep = "")

