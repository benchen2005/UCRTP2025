###############################################################################
# GO bar & dot plots — BIG LEFT LABELS (Arial) + distinct colors per ontology
###############################################################################

# ==== USER SETTINGS ====
tsv_file    <- "/Users/ben/Documents/UCRTP/Filtered Events/Filtered SE.tsv"
psi_cutoff  <- 0.10
fdr_cutoff  <- 0.05
out_prefix  <- "/Users/ben/Documents/UCRTP/GO/SE_GO"   # folder+basename

show_n      <- 15     # top N terms to show
base_font   <- 24     # general UI text size
label_font  <- 44     # *** LEFT y-axis labels (GO terms) ***
value_font  <- 14     # numbers on bars (keep smaller so labels stand out)
wrap_width  <- 60     # wrap GO names to this many chars per line

# ==== PACKAGES ====
pkgs <- c("dplyr","readr","clusterProfiler","org.Hs.eg.db",
          "ggplot2","scales","stringr")
need <- setdiff(pkgs, rownames(installed.packages()))xx
if (length(need)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(need, update = FALSE, ask = FALSE)
}
invisible(lapply(pkgs, library, character.only = TRUE))

# ==== READ & FILTER ====
dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

numify <- function(x) if (is.character(x)) readr::parse_number(x) else as.numeric(x)

se <- read.delim(tsv_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("geneSymbol","PSI","FDR") %in% names(se)))

se <- se |>
  dplyr::mutate(PSI = numify(PSI), FDR = numify(FDR))

# avoid Inf in -log10(FDR)
if (any(se$FDR == 0, na.rm = TRUE)) {
  eps <- min(se$FDR[se$FDR > 0], na.rm = TRUE) * 0.5
  se$FDR[se$FDR == 0] <- eps
}

sig_genes <- se |>
  dplyr::filter(is.finite(PSI), is.finite(FDR),
                abs(PSI) > psi_cutoff, FDR < fdr_cutoff) |>
  dplyr::pull(geneSymbol) |>
  unique()

cat("Significant genes:", length(sig_genes), "\n")
if (!length(sig_genes)) stop("No genes passed thresholds; adjust psi_cutoff/fdr_cutoff.")

# ==== ID mapping ====
entrez <- clusterProfiler::bitr(sig_genes,
                                fromType = "SYMBOL",
                                toType   = "ENTREZID",
                                OrgDb    = org.Hs.eg.db)
entrez_ids <- unique(entrez$ENTREZID)
cat("Mapped to Entrez IDs:", length(entrez_ids), "\n")
if (length(entrez_ids) < 5) stop("Too few mapped Entrez IDs (n < 5).")

# ==== GO enrichment ====
ego_args <- list(OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                 pAdjustMethod = "BH", pvalueCutoff = 1)

go_bp <- do.call(enrichGO, c(list(gene = entrez_ids, ont = "BP"), ego_args))
go_cc <- do.call(enrichGO, c(list(gene = entrez_ids, ont = "CC"), ego_args))
go_mf <- do.call(enrichGO, c(list(gene = entrez_ids, ont = "MF"), ego_args))

# save full tables
write.csv(as.data.frame(go_bp), paste0(out_prefix, "_BP.csv"), row.names = FALSE)
write.csv(as.data.frame(go_cc), paste0(out_prefix, "_CC.csv"), row.names = FALSE)
write.csv(as.data.frame(go_mf), paste0(out_prefix, "_MF.csv"), row.names = FALSE)

# ==== Prep & Theme helpers ====
prep_go_df <- function(go_obj, top_n = 15, wrap_width = 60) {
  df <- as.data.frame(go_obj)
  if (is.null(df) || !nrow(df)) return(NULL)
  
  # GeneRatio "12/153" -> numeric
  df$GeneRatioNum <- sapply(strsplit(df$GeneRatio, "/"),
                            function(v) as.numeric(v[1]) / as.numeric(v[2]))
  df$neglog10FDR <- -log10(pmax(df$p.adjust, .Machine$double.xmin))
  
  # order by FDR then Count; take top_n
  df <- df[order(df$p.adjust, -df$Count), , drop = FALSE]
  n_show <- min(top_n, nrow(df))
  df <- df[seq_len(n_show), , drop = FALSE]
  
  # Wrapped label + counts to aid readability
  df$Label <- paste0(stringr::str_wrap(df$Description, width = wrap_width),
                     "  (", df$Count, ")")
  df$Label <- factor(df$Label, levels = rev(df$Label))
  df
}

theme_big_arial <- function() {
  ggplot2::theme_minimal(base_size = base_font, base_family = "Arial") +
    ggplot2::theme(
      text          = ggplot2::element_text(family = "Arial"),
      plot.title    = ggplot2::element_text(size = base_font + 8, face = "bold"),
      axis.title.x  = ggplot2::element_text(size = base_font + 4),
      axis.text.x   = ggplot2::element_text(size = base_font - 2),
      axis.text.y   = ggplot2::element_text(size = label_font, lineheight = 0.95), # <<<< BIG LEFT LABELS
      legend.title  = ggplot2::element_text(size = base_font),
      legend.text   = ggplot2::element_text(size = base_font - 2),
      plot.margin   = ggplot2::margin(20, 40, 20, 100)  # extra left margin so big labels don't clip
    )
}

# Distinct colors per ontology
bar_fill <- function(tag) switch(tag,
                                 BP = "#1f78b4",   # blue
                                 CC = "#33a02c",   # green
                                 MF = "#e31a1c",   # red
                                 "#6a3d9a"         # fallback
)
dot_grad <- function(tag) switch(tag,
                                 BP = c(low = "#c6dbef", high = "#08306b"),
                                 CC = c(low = "#c7e9c0", high = "#005a32"),
                                 MF = c(low = "#fdd0a2", high = "#a63603"),
                                 c(low = "#dadaeb", high = "#3f007d")
)

# ==== Plotters (BIG left labels, modest right numbers) ====
plot_go_bar <- function(df, tag, title_txt, file_png,
                        w = 20, h = 0.9*show_n + 5) {
  if (is.null(df)) return(invisible())
  g <- ggplot2::ggplot(df, ggplot2::aes(x = neglog10FDR, y = Label)) +
    ggplot2::geom_col(fill = bar_fill(tag), width = 0.72) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", neglog10FDR)),
                       hjust = -0.12, size = value_font, family = "Arial") +   # keep small
    ggplot2::coord_cartesian(xlim = c(0, max(df$neglog10FDR) * 1.25)) +
    ggplot2::labs(title = title_txt,
                  x = expression(-log[10]~"FDR (BH)"),
                  y = NULL) +
    theme_big_arial()
  ggplot2::ggsave(file_png, g, width = w, height = h, dpi = 300)
}

plot_go_dot <- function(df, tag, title_txt, file_png,
                        w = 20, h = 0.9*show_n + 5) {
  if (is.null(df)) return(invisible())
  grad <- dot_grad(tag)
  g <- ggplot2::ggplot(df, ggplot2::aes(x = neglog10FDR, y = Label)) +
    ggplot2::geom_point(ggplot2::aes(size = Count, color = GeneRatioNum)) +
    ggplot2::scale_color_gradient(low = grad["low"], high = grad["high"],
                                  name = "Gene ratio") +
    ggplot2::scale_size(range = c(8, 20), name = "Hits") +
    ggplot2::labs(title = title_txt,
                  x = expression(-log[10]~"FDR (BH)"),
                  y = NULL) +
    theme_big_arial()
  ggplot2::ggsave(file_png, g, width = w, height = h, dpi = 300)
}

# ==== Build & save all plots ====
bp_df <- prep_go_df(go_bp, show_n, wrap_width)
cc_df <- prep_go_df(go_cc, show_n, wrap_width)
mf_df <- prep_go_df(go_mf, show_n, wrap_width)

plot_go_bar(bp_df, "BP", "GO Biological Process (top terms)",
            paste0(out_prefix, "_BP_bar.png"))
plot_go_dot(bp_df, "BP", "GO Biological Process (top terms)",
            paste0(out_prefix, "_BP_dot.png"))

plot_go_bar(cc_df, "CC", "GO Cellular Component (top terms)",
            paste0(out_prefix, "_CC_bar.png"))
plot_go_dot(cc_df, "CC", "GO Cellular Component (top terms)",
            paste0(out_prefix, "_CC_dot.png"))

plot_go_bar(mf_df, "MF", "GO Molecular Function (top terms)",
            paste0(out_prefix, "_MF_bar.png"))
plot_go_dot(mf_df, "MF", "GO Molecular Function (top terms)",
            paste0(out_prefix, "_MF_dot.png"))

cat("\n✔ Saved figures with HUGE left labels (Arial) & distinct colors:\n  ",
    normalizePath(paste0(out_prefix, "_BP_bar.png")), "\n  ",
    normalizePath(paste0(out_prefix, "_BP_dot.png")), "\n  ",
    normalizePath(paste0(out_prefix, "_CC_bar.png")), "\n  ",
    normalizePath(paste0(out_prefix, "_CC_dot.png")), "\n  ",
    normalizePath(paste0(out_prefix, "_MF_bar.png")), "\n  ",
    normalizePath(paste0(out_prefix, "_MF_dot.png")), "\n", sep = "")

