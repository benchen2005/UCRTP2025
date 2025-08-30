library(ggplot2)
setwd("/Users/ben/Downloads/Filtered SE.tsv")
f <- "Filtered SE - SE_ FDR _ 0.001, _PSI_ _ 0.1.tsv"                  # adjust path if needed
stopifnot(file.exists(f))                   # clear error if file missing

se <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Coerce just in case
se$IncLevelDifference <- as.numeric(se$IncLevelDifference)
se$FDR                <- as.numeric(se$FDR)

volcano <- ggplot(se, aes(x = IncLevelDifference,
                          y = -log10(FDR))) +
  geom_point(aes(colour = abs(IncLevelDifference) > 0.1 & FDR < 0.001),
             alpha = 0.7) +
  scale_colour_manual(values = c("grey50", "red"),
                      labels  = c("NS", "Sig"),
                      name    = "ΔPSI > 0.1 & FDR < 0.001") +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(x = "PSI difference (ΔPSI)",
       y = expression(-log[10]~FDR),
       title = "Volcano plot – Skipped-exon events") +
  theme_minimal()

print(volcano)




# Optional: save to file
ggplot2::ggsave("volcano_SE.png", volcano,
                width = 8, height = 6, dpi = 300)


##### 1.  Load the SE table (already in your environment) ######################
se <- read.table("SE.MATS.JCEC.DN.txt",
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE)

##### 2.  Filter for “hits” ####################################################
library(dplyr)   # if not attached yet

hits <- se %>%
  filter(abs(IncLevelDifference) > 0.10,      # effect-size cut-off
         FDR < 0.05)                          # significance cut-off

##### 3.  Inspect & count ######################################################
nrow(hits)           # how many significant SE events?
head(hits, 3)        # peek at first few rows

##### 4.  Save to file (tab-separated) #########################################
write.table(hits,
            file = "SE_hits_DeltaPSI0.1_FDR0.05.txt",
            sep  = "\t",
            quote = FALSE,
            row.names = FALSE)

top20 <- hits %>% arrange(desc(abs(IncLevelDifference))) %>% slice(1:20)

