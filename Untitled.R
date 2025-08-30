# ---- Inputs ----
len <- 580
features <- data.frame(
  start = c(19, 350, 490),
  end   = c(124, 380, 537),
  name  = c("XTBD", "DBD", "RBD")
)

library(ggplot2)
library(dplyr)

# Ticks at boundaries + termini
ticks <- tibble::tibble(x = c(0, features$start, features$end, len)) %>% distinct()

# ---- Plot ----
p <- ggplot() +
  # main protein bar
  geom_rect(aes(xmin = 0, xmax = len, ymin = 0.45, ymax = 0.55),
            fill = "#F6B27B", color = NA) +
  # domain separators (thin)
  geom_segment(data = ticks, aes(x = x, xend = x, y = 0.42, yend = 0.58),
               linewidth = 0.35, color = "#3A3A3A") +
  # domain labels above
  geom_text(data = features,
            aes(x = (start + end)/2, y = 0.70, label = name),
            fontface = "bold", size = 4, color = "#2F2F2F") +
  # residue ranges under each domain (optional)
  geom_text(data = features,
            aes(x = (start + end)/2, y = 0.28,
                label = paste0("(", start, "–", end, " aa)")),
            size = 3.2, color = "#4A4A4A") +
  # title (gene & alias)
  annotate("text", x = 0, y = 0.85, hjust = 0,
           label = "CDKN2AIP (CARF)", fontface = "bold", size = 5.2) +
  # termini labels
  annotate("text", x = 0,   y = 0.15, label = "N (1 aa)",  size = 3.2, hjust = 0) +
  annotate("text", x = len, y = 0.15, label = "C (580 aa)", size = 3.2, hjust = 1) +
  # clean theme
  coord_cartesian(xlim = c(0, len), ylim = c(0, 1)) +
  theme_void() +
  theme(plot.margin = margin(20, 28, 20, 28))

# Preview and export
print(p)
ggsave("CDKN2AIP_domains.svg", p, width = 8, height = 2.2)   # vector
ggsave("CDKN2AIP_domains.pdf", p, width = 8, height = 2.2)
ggsave("CDKN2AIP_domains.png", p, width = 8, height = 2.2, dpi = 600)

