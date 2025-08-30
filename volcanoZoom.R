###############################################################################
# 0 ── USER SETTINGS                                                          #
###############################################################################
file_path   <- "/Users/ben/Documents/UCRTP/Filtered Events/Filtered SE.tsv"
psi_thresh  <- 0.10                    # |ΔPSI| threshold
fdr_thresh  <- 0.001                   # FDR threshold
genes_interest <- c(
  "SRSF5","SRSF2","SRSF4","SRSF6","SRSF3","SRSF11","BRCA1","RSRP1","Sirt2",
  "Sirt6","NAMPT","HNRNPA1","BCL2L1","BCL2","SLC7A11","UBE2F","ESRP1",
  "HNRNPD","HNRNPU","UBE3A","UBE4B","NOTCH1","SIRT7","HRAS","PFKL","XBP1"
)

###############################################################################
# 1 ── LOAD PACKAGES (auto-install)                                           #
###############################################################################
pkgs <- c("ggplot2","ggrepel","readr","dplyr","patchwork")
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need)
lapply(pkgs, library, character.only = TRUE)

###############################################################################
# 2 ── READ & TIDY                                                            #
###############################################################################
stopifnot(file.exists(file_path))
se <- read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)

# ensure required columns
stopifnot(all(c("geneSymbol","PSI","FDR") %in% names(se)))

numify <- function(x) if (is.character(x)) readr::parse_number(x) else as.numeric(x)
se <- se %>% mutate(PSI = numify(PSI), FDR = numify(FDR))

###############################################################################
# 3 ── LABEL & PREP                                                           #
###############################################################################
se <- se %>%
  mutate(Sig   = abs(PSI) > psi_thresh & FDR < fdr_thresh,
         GOI   = geneSymbol %in% genes_interest,
         label = ifelse(GOI, geneSymbol, NA))

# handle FDR == 0 so log10 is finite
min_non0 <- min(se$FDR[se$FDR > 0], na.rm = TRUE)
se$FDR[se$FDR == 0] <- min_non0 * 0.5

y_cap <- quantile(-log10(se$FDR), 0.99, na.rm = TRUE) * 1.10  # full-plot ceiling
top_band <- quantile(-log10(se$FDR), 0.95, na.rm = TRUE)      # zoom start

###############################################################################
# 4 ── BASE VOLCANO PLOT                                                      #
###############################################################################
base_plot <- ggplot(se, aes(PSI, -log10(FDR))) +
  geom_point(data = subset(se, !Sig), colour = "grey75",  size = 1) +
  geom_point(data = subset(se,  Sig & !GOI), colour = "#F8766D", size = 1.2) +
  geom_point(data = subset(se,  GOI), shape = 21, stroke = .4,
             fill = "dodgerblue2", colour = "dodgerblue3", size = 3) +
  geom_text_repel(data = subset(se, GOI), aes(label = label), fontface = "bold",
                  box.padding = .3, point.padding = .15, max.overlaps = Inf,
                  min.segment.length = 0, size = 3) +
  geom_vline(xintercept = c(-psi_thresh, psi_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_thresh),          linetype = "dashed") +
  scale_y_continuous(limits = c(0, y_cap),
                     expand = expansion(mult = c(0, .02))) +
  labs(title = "Volcano plot – Skipped-exon (SE) events",
       subtitle = paste("Red = Sig (|ΔPSI|>", psi_thresh,
                        "& FDR<", fdr_thresh, ")   •   Blue = Genes of interest"),
       x = "PSI difference (ΔPSI)",
       y = expression(-log[10]~FDR)) +
  theme_minimal(base_size = 11)

###############################################################################
# 5 ── ZOOM PANEL (top 5 %)                                                   #
###############################################################################
zoom_plot <- base_plot +
  coord_cartesian(ylim = c(top_band, y_cap)) +
  labs(subtitle = paste0("Zoom: top 5 %  (–log10 FDR ≥ ",
                         round(top_band,1), ")")) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())

###############################################################################
# 6 ── COMBINE & SAVE                                                         #
###############################################################################
combo <- base_plot / zoom_plot + plot_layout(heights = c(3, 1))

ggsave("volcano_SE_GOI_combo.png", combo,
       width = 9, height = 10, dpi = 300)
ggsave("volcano_SE_GOI_zoom.png",  zoom_plot,
       width = 7, height = 5, dpi = 300)

print(combo)

cat("\nSaved:\n  • volcano_SE_GOI_combo.png  (full + zoom)\n",
    " • volcano_SE_GOI_zoom.png   (zoom only)\n")

