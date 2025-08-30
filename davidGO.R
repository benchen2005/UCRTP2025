install.packages("showtext")        # if not yet installed
library(showtext)
font_add("Arial", regular = "Arial.ttf")  # system Arial
showtext_auto()                  # render every plot with showtext

###############################################################################
# 0 ── USER-EDITABLE PATHS                                                    #
###############################################################################
root_dir   <- "/Users/ben/Documents/UCRTP/Filtered Events"
out_root   <- "/Users/ben/Documents/UCRTP/GO_all_events"
event_files <- c(
  SE   = "noFilterFDR.SE.tsv",
  RI   = "noFilterFDR.RI.tsv",
  A3SS = "noFilterFDR.A3SS.tsv",
  A5SS = "noFilterFDR.A5SS.tsv",
  MXE  = "noFilterFDR.MXE.tsv"
)

psi_cutoff <- 0.10
fdr_cutoff <- 0.01
###############################################################################

###############################################################################
# 1 ── PACKAGES & FONT (Arial via showtext)                                   #
###############################################################################
pkgs <- c("dplyr","readr","clusterProfiler","org.Hs.eg.db",
          "enrichplot","ggplot2","showtext")
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(need, ask = FALSE, update = FALSE)
}
lapply(pkgs, library, character.only = TRUE)

## --- load Arial once -------------------------------------------------------
font_add("Arial", regular = "Arial.ttf")
showtext_auto()

###############################################################################
# 2 ── BAR-PLOT FUNCTION (DAVID-style)                                        #
###############################################################################
plot_bar <- function(go_obj, n = 10, title = "", bar_colour = "#1f78b4",
                     file) {
  
  df <- as.data.frame(go_obj) %>%
    slice_head(n = n) %>%
    mutate(Description = paste0(Description, " (", Count, ")")) %>%
    arrange(p.adjust)
  df$Description <- factor(df$Description, levels = rev(df$Description))
  
  g <- ggplot(df, aes(x = -log10(p.adjust), y = Description)) +
    geom_col(fill = bar_colour, width = 0.7) +
    scale_x_continuous(trans = "log10",
                       breaks = c(1e-1,1e-2,1e-3,1e-4,1e-6,1e-9),
                       labels = scales::trans_format(
                         "log10", scales::math_format(10^.x))) +
    labs(title = title,
         x     = expression(italic("p") * "-value"),
         y     = NULL) +
    theme_minimal(base_size = 14, base_family = "Arial") +
    theme(
      axis.title.x = element_text(size = 14, margin = margin(t = 8)),
      axis.text.y  = element_text(size = 12, family = "Arial"),
      axis.text.x  = element_text(size = 12, family = "Arial"),
      plot.title   = element_text(size = 16, face = "bold", family = "Arial"),
      plot.margin  = margin(10, 40, 15, 10)
    )
  
  ggsave(file, g,
         width  = 8.5,
         height = 0.4 * n + 2.7,   # extra space for axis title
         dpi    = 300)
}

###############################################################################
# 3 ── LOOP OVER EVENT TYPES                                                  #
###############################################################################
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

for (etype in names(event_files)) {
  
  message("\n=====  Processing ", etype, "  =====")
  tsv_path <- file.path(root_dir, event_files[[etype]])
  stopifnot(file.exists(tsv_path))
  
  out_dir <- file.path(out_root, etype)
  dir.create(out_dir, showWarnings = FALSE)
  
  # --- read & numeric -------------------------------------------------------
  tab <- read.delim(tsv_path, check.names = FALSE, stringsAsFactors = FALSE)
  numify <- function(x) if (is.character(x)) readr::parse_number(x) else as.numeric(x)
  tab <- tab %>% mutate(PSI = numify(PSI), FDR = numify(FDR))
  
  # --- significant events ---------------------------------------------------
  sig <- tab %>% filter(abs(PSI) > psi_cutoff, FDR < fdr_cutoff)
  
  write.csv(sig, file.path(out_dir, paste0(etype, "_sig_events.csv")), row.names = FALSE)
  
  genes <- unique(sig$geneSymbol)
  writeLines(genes, file.path(out_dir, paste0(etype, "_sig_genes.txt")))
  cat("  Significant events:", nrow(sig), "  genes:", length(genes), "\n")
  
  if (length(genes) < 2) { cat("  Too few genes — skipped.\n"); next }
  
  # --- SYMBOL → ENTREZ ------------------------------------------------------
  entrez <- suppressMessages(
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db))
  entrez_ids <- unique(entrez$ENTREZID)
  if (length(entrez_ids) < 2) { cat("  <2 Entrez IDs — skipped.\n"); next }
  
  # --- GO enrichment --------------------------------------------------------
  go_bp <- enrichGO(entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP",
                    pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE)
  go_cc <- enrichGO(entrez_ids, OrgDb = org.Hs.eg.db, ont = "CC",
                    pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE)
  go_mf <- enrichGO(entrez_ids, OrgDb = org.Hs.eg.db, ont = "MF",
                    pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE)
  
  write.csv(go_bp, file.path(out_dir, paste0(etype, "_GO_BP.csv")), row.names = FALSE)
  write.csv(go_cc, file.path(out_dir, paste0(etype, "_GO_CC.csv")), row.names = FALSE)
  write.csv(go_mf, file.path(out_dir, paste0(etype, "_GO_MF.csv")), row.names = FALSE)
  
  # --- bar plots ------------------------------------------------------------
  plot_bar(go_cc, 8,  paste0("DAVID GOTERM – CC (", etype, ")"),
           bar_colour = "#1b9e77",
           file = file.path(out_dir, paste0(etype, "_GO_CC_bar.png")))
  plot_bar(go_mf, 8,  paste0("DAVID GOTERM – MF (", etype, ")"),
           bar_colour = "#d95f02",
           file = file.path(out_dir, paste0(etype, "_GO_MF_bar.png")))
  plot_bar(go_bp, 10, paste0("DAVID GOTERM – BP (", etype, ")"),
           bar_colour = "#1f78b4",
           file = file.path(out_dir, paste0(etype, "_GO_BP_bar.png")))
  
  cat("  ✔ Outputs written to ", normalizePath(out_dir), "\n")
}

cat("\nAll event types finished. Results root: ", normalizePath(out_root), "\n")
