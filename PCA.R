###############################################################################
# PCA of samples using PSI values from rMATS table
# - Reads your TSV
# - Builds PSI matrix (events × samples)
# - Filters low-information rows
# - Runs PCA (samples are observations)
# - Saves PC1–PC2 scatter + scree plot
###############################################################################

# ===== USER SETTINGS =====
file_path <- "/Users/ben/Documents/UCRTP/Filtered Events/noFilterFDR.SE.tsv"
psi_cols  <- c("Ctrl1","Ctrl2","Ctrl3","KD1","KD2","KD3")  # PSI columns that exist in your TSV
out_dir   <- "plots"
out_png   <- file.path(out_dir, "PCA_PSI_PC1_PC2.png")
out_pdf   <- file.path(out_dir, "PCA_PSI_PC1_PC2.pdf")
out_scree <- file.path(out_dir, "PCA_scree.png")

# filtering knobs (optional, good defaults)
min_non_na_per_row <- 4     # drop events with fewer than this many non-NA PSI values
top_n_by_var       <- 3000  # keep top N most variable events (set Inf for all)

# ===== PACKAGES (namespaced to avoid conflicts) =====
pkgs <- c("dplyr","readr","ggplot2","ggrepel","scales")
new  <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new)
invisible(lapply(pkgs, library, character.only = TRUE))

numify <- function(x) if (is.character(x)) readr::parse_number(x) else as.numeric(x)

# ===== LOAD TSV & CHECK COLUMNS =====
stopifnot(file.exists(file_path))
tab <- utils::read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)

missing <- setdiff(psi_cols, names(tab))
if (length(missing)) {
  stop("Missing PSI columns: ", paste(missing, collapse = ", "),
       "\nPresent columns:\n", paste(names(tab), collapse = ", "))
}

# ===== BUILD PSI MATRIX (events × samples) =====
psi_df <- tab |>
  dplyr::mutate(event_id = paste0(
    if ("geneSymbol" %in% names(tab)) tab$geneSymbol else "event",
    "_", dplyr::row_number()
  )) |>
  dplyr::select(event_id, dplyr::all_of(psi_cols))

# unique rownames
rownames(psi_df) <- psi_df$event_id
psi_df$event_id  <- NULL

# coerce to numeric
psi_df[] <- lapply(psi_df, numify)

# drop rows with too many NAs
keep <- rowSums(is.na(psi_df)) <= (ncol(psi_df) - min_non_na_per_row)
psi_df <- psi_df[keep, , drop = FALSE]

# if desired, keep top-N by variance (helps stabilize PCA)
if (is.finite(top_n_by_var) && nrow(psi_df) > top_n_by_var) {
  vars <- apply(psi_df, 1, stats::var, na.rm = TRUE)
  psi_df <- psi_df[order(vars, decreasing = TRUE)[1:top_n_by_var], , drop = FALSE]
}

# impute any remaining NA by row mean (simple & safe for PCA)
row_means <- rowMeans(psi_df, na.rm = TRUE)
na_idx <- which(is.na(psi_df), arr.ind = TRUE)
if (nrow(na_idx)) {
  psi_df[na_idx] <- row_means[na_idx[, "row"]]
}

# ===== PCA (samples are observations) =====
# prcomp expects rows = observations (samples), cols = variables (events)
X <- t(as.matrix(psi_df))        # now samples × events
stopifnot(all(is.finite(X)))

pc <- stats::prcomp(X, center = TRUE, scale. = TRUE)
expl_var <- pc$sdev^2
expl_var <- expl_var / sum(expl_var)

pc_scores <- as.data.frame(pc$x)  # samples × PCs
pc_scores$Sample <- rownames(pc_scores)

# sample groups from names: Ctrl* vs KD*
pc_scores$Group <- ifelse(grepl("^Ctrl", pc_scores$Sample), "Control", "KnockDown")
pc_scores$Group <- factor(pc_scores$Group, levels = c("Control","KnockDown"))

pc1_lab <- paste0("PC1 (", scales::percent(expl_var[1], accuracy = 0.1), ")")
pc2_lab <- paste0("PC2 (", scales::percent(expl_var[2], accuracy = 0.1), ")")

message("PCA built from ", nrow(psi_df), " events and ", ncol(psi_df), " samples.")

# ===== PLOTS =====
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# PC1 vs PC2 scatter with ellipses + labels
p <- ggplot2::ggplot(pc_scores, ggplot2::aes(PC1, PC2, color = Group, shape = Group)) +
  ggplot2::geom_point(size = 4, alpha = 0.9) +
  ggplot2::stat_ellipse(type = "norm", level = 0.68, linewidth = 0.7, alpha = 0.25) +
  ggrepel::geom_text_repel(ggplot2::aes(label = Sample),
                           size = 5.2, fontface = "bold", seed = 42,
                           box.padding = 0.35, point.padding = 0.25, min.segment.length = 0) +
  ggplot2::scale_color_manual(values = c(Control = "#2b8cbe", KnockDown = "#d95f02")) +
  ggplot2::labs(title = "PCA of PSI profiles",
                subtitle = paste0("Built from ", nrow(psi_df), " SE events (filtered & scaled)"),
                x = pc1_lab, y = pc2_lab, color = "Group", shape = "Group") +
  ggplot2::theme_minimal(base_size = 16) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 18),
                 panel.grid.minor = ggplot2::element_blank())

ggplot2::ggsave(out_png, p, width = 9.5, height = 7.5, dpi = 300)
ggplot2::ggsave(out_pdf, p, width = 9.5, height = 7.5)
message("✔ PCA PC1–PC2 saved to: ", normalizePath(out_png))

# Scree plot (variance explained)
scree_df <- data.frame(PC = paste0("PC", seq_along(expl_var)),
                       var = expl_var,
                       cum = cumsum(expl_var))

s <- ggplot2::ggplot(scree_df[1:10, ], ggplot2::aes(x = PC, y = var)) +
  ggplot2::geom_col(fill = "grey40") +
  ggplot2::geom_line(ggplot2::aes(y = cum, group = 1), linewidth = 1.1, color = "#2b8cbe") +
  ggplot2::geom_point(ggplot2::aes(y = cum), size = 2.5, color = "#2b8cbe") +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::labs(title = "PCA scree plot (top 10 PCs)",
                x = NULL, y = "Variance explained") +
  ggplot2::theme_minimal(base_size = 16) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 18))
ggplot2::ggsave(out_scree, s, width = 8.5, height = 6.5, dpi = 300)
message("✔ Scree plot saved to: ", normalizePath(out_scree))
