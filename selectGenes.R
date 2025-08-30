###############################################################################
# 0 ── user settings                                                          #
###############################################################################
file_path   <- "/Users/ben/Documents/Filtered SE.tsv"   #  <-- adjust if needed
psi_cutoff  <- 0.80                                     #  |ΔPSI| threshold
out_file    <- "genes_absPSI_gt_0.5.txt"                #  output list

###############################################################################
# 1 ── install / load packages                                                #
###############################################################################
pkgs <- c("readr", "dplyr")
new  <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new)) install.packages(new)
invisible(lapply(pkgs, library, character.only = TRUE))

###############################################################################
# 2 ── read & clean the TSV                                                   #
###############################################################################
stopifnot(file.exists(file_path))

se <- read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)

req_cols <- c("geneSymbol", "PSI")  # FDR not needed for this task
if (!all(req_cols %in% names(se))) {
  stop("The file must contain columns: ", paste(req_cols, collapse = ", "))
}

# robust numeric conversion
clean_num <- function(x) {
  if (is.character(x)) readr::parse_number(x) else as.numeric(x)
}

se$PSI <- clean_num(se$PSI)

###############################################################################
# 3 ── extract genes with |ΔPSI| > 0.75                                        #
###############################################################################
genes_big_shift <- se %>%
  filter(abs(PSI) > psi_cutoff) %>%   # at least one event meets the cutoff
  distinct(geneSymbol) %>%            # unique genes
  arrange(geneSymbol) %>%             # alphabetical order
  pull(geneSymbol)

###############################################################################
# 4 ── output                                                                 #
###############################################################################
# write to disk
writeLines(genes_big_shift, out_file)

# print summaries
cat("\nGenes with at least one SE event where |ΔPSI| >", psi_cutoff, ":\n",
    paste(genes_big_shift, collapse = ", "), "\n\n",
    "Total genes:", length(genes_big_shift), "\n",
    "Saved to:", normalizePath(out_file), "\n")
