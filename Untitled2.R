###############################################################################
# 0 ── user settings                                                          #
###############################################################################
file_path   <- "/Users/ben/Documents/Filtered SE.tsv"   # ← TSV file
top_n       <- 50        # top N most-variable events (Inf for all)
ctrl_regex  <- "^Ctrl"   # pattern for Control sample columns
kd_regex    <- "^KD"     # pattern for Knock-down sample columns
out_png     <- "SE_PSI_heatmap.png"

###############################################################################
# 1 ── packages (auto-install if missing)                                     #
###############################################################################
pkgs <- c("dplyr","tidyr","readr","pheatmap")
new  <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new)
invisible(lapply(pkgs, library, character.only = TRUE))

###############################################################################
# 2 ── read TSV & detect per-sample PSI columns                               #
###############################################################################
stopifnot(file.exists(file_path))
se <- read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)

psi_cols <- names(se)[
  sapply(se, \(col) {
    is_char <- is.character(col) || is.factor(col)
    looks_numeric <- all(grepl("^[0-9., ]*$", col))
    has_comma_or_dot <- any(grepl("[.,]", col))
    is_char && looks_numeric && has_comma_or_dot
  })
]

if (length(psi_cols) < 2)
  stop("Detected fewer than two PSI columns.\n",
       "Please edit `psi_cols` manually after running `names(se)`.")

message("PSI columns detected: ", paste(psi_cols, collapse = ", "))

###############################################################################
# 3 ── widen PSI matrix (one column per replicate)                            #
###############################################################################
psi_long <- se %>%
  mutate(event_id = row_number()) %>%
  mutate(across(all_of(psi_cols), as.character)) %>%    # ensure single type
  select(event_id, all_of(psi_cols)) %>%
  pivot_longer(-event_id, names_to = "Group", values_to = "PSI_str") %>%
  separate_rows(PSI_str, sep = ",") %>%                 # split replicates
  group_by(event_id, Group) %>%
  mutate(rep = row_number()) %>% ungroup() %>%
  mutate(Sample = paste0(Group, "_", rep)) %>%
  select(event_id, Sample, PSI_str)

psi_wide <- psi_long %>%
  pivot_wider(names_from = Sample, values_from = PSI_str) %>%
  column_to_rownames("event_id") %>%
  mutate(across(everything(), readr::parse_number))

psi_mat <- as.matrix(psi_wide)

###############################################################################
# 4 ── keep top-N most-variable events                                        #
###############################################################################
if (is.finite(top_n) && nrow(psi_mat) > top_n) {
  vars <- apply(psi_mat, 1, var, na.rm = TRUE)
  psi_mat <- psi_mat[order(vars, decreasing = TRUE)[1:top_n], ]
}

###############################################################################
# 5 ── build sample annotation                                                #
###############################################################################
sample_groups <- ifelse(grepl(ctrl_regex, colnames(psi_mat)), "Control",
                        ifelse(grepl(kd_regex, colnames(psi_mat)), "KnockDown",
                               "Other"))
anno <- data.frame(Group = factor(sample_groups, levels = unique(sample_groups)))
rownames(anno) <- colnames(psi_mat)

###############################################################################
# 6 ── draw heat-map                                                          #
###############################################################################
hm <- pheatmap(
  psi_mat,
  annotation_col = anno,
  show_rownames  = FALSE,
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  color          = colorRampPalette(c("royalblue","white","firebrick3"))(100),
  main           = paste0("PSI heat-map (", nrow(psi_mat),
                          " events; top ", ifelse(is.finite(top_n), top_n, "all"), ")")
)

###############################################################################
# 7 ── save                                                                   #
###############################################################################
ggsave(out_png, hm, width = 10, height = 8, dpi = 300)
cat("Heat-map saved to:", normalizePath(out_png), "\n")
