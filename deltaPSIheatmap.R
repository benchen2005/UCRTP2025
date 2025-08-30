###############################################################################
# USER SETTINGS                                                               #
###############################################################################
file_path <- "/Users/ben/Documents/Filtered SE.tsv"              # TSV path
ctrl_cols <- c("Ctrl1","Ctrl2","Ctrl3")                          # control cols
kd_cols   <- c("KD1","KD2","KD3")                                # KD cols
out_dir   <- "plots"
out_png   <- file.path(out_dir, "SE_DeltaPSI_pairwise.png")

###############################################################################
# PACKAGES                                                                    #
###############################################################################
pkgs <- c("dplyr","pheatmap")
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need)
lapply(pkgs, library, character.only = TRUE)

###############################################################################
# READ TSV & COERCE SAMPLE COLUMNS                                            #
###############################################################################
stopifnot(file.exists(file_path))
se <- read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)

missing <- setdiff(c(ctrl_cols, kd_cols), names(se))
if (length(missing))
  stop("Missing columns: ", paste(missing, collapse = ", "))

se <- se %>%
  mutate(across(all_of(c(ctrl_cols, kd_cols)),
                ~ as.numeric(as.character(.))))  # ensure numeric

###############################################################################
# COMPUTE ΔΨ (KD – Control) FOR EACH PAIR                                     #
###############################################################################
delta_df <- se %>%
  mutate(event_id = paste0(geneSymbol, "_", row_number()),
         `ΔΨ KD1−Ctrl1` = .data[[kd_cols[1]]] - .data[[ctrl_cols[1]]],
         `ΔΨ KD2−Ctrl2` = .data[[kd_cols[2]]] - .data[[ctrl_cols[2]]],
         `ΔΨ KD3−Ctrl3` = .data[[kd_cols[3]]] - .data[[ctrl_cols[3]]]) %>%
  select(event_id, starts_with("ΔΨ"))

rownames(delta_df) <- delta_df$event_id
delta_df$event_id  <- NULL
delta_mat <- as.matrix(delta_df)

###############################################################################
# COLOUR SCALE (−1 to +1)                                                     #
###############################################################################
breaks <- seq(-1, 1, length.out = 101)
cols   <- colorRampPalette(c("royalblue", "white", "firebrick3"))(100)

###############################################################################
# OUTPUT FOLDER                                                               #
###############################################################################
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# DRAW HEAT-MAP AND SAVE PNG                                                  #
###############################################################################
title_main <- "Pair-wise ΔΨ heat-map  (KD – Control)"
title_sub  <- "Firebrick = More included in KD · Royal-blue = Less included in KD"

png(out_png, width = 2400, height = 1800, res = 300)
pheatmap::pheatmap(
  delta_mat,
  cluster_rows   = TRUE,
  cluster_cols   = FALSE,
  show_rownames  = FALSE,
  color          = cols,
  breaks         = breaks,
  main           = paste0(title_main, "\n", title_sub),
  angle_col      = 0                                # keep column labels horizontal
)
dev.off()

# Show interactively as well
pheatmap::pheatmap(
  delta_mat,
  cluster_rows  = TRUE,
  cluster_cols  = FALSE,
  show_rownames = FALSE,
  color         = cols,
  breaks        = breaks,
  main          = paste0(title_main, "\n", title_sub),
  angle_col     = 0
)

cat("✔ Heat-map saved to", normalizePath(out_png), "\n")
