library(tidyverse)

dat <- read_csv("raw_data_for_lineage.csv")

library(dplyr)
library(tidyr)
library(ggplot2)

# Make Cluster a factor so boxes are discrete
dat <- dat %>%
  mutate(Cluster = factor(Cluster, levels = 1:12))

mutation_cols <- names(dat)[-(1:2)]  # or explicitly c("D3L","P13L",...)

# Reshape all escape columns to long format
dat_long <- dat %>%
  pivot_longer(
    cols = -c(Antibody, Cluster),   # all escape columns
    names_to = "Mutation",
    values_to = "Escape"
  ) %>%
  mutate(
    Mutation = factor(Mutation, levels = mutation_cols)  # preserve order
  )

custom_colors_14 <- c("#4DAF4A", "#E4C899", "#377EB8", "#E41A1C", 
                      "#A65628", "#984EA3", "grey40", "#D3AF37",
                      # Extensions (complementary tones)
                      "#F781BF", "#66C2A5", "#999999", 
                      "#FF7F00")

# Boxplots: one panel per mutation, boxes = clusters
ggplot(dat_long, aes(x = Cluster, y = Escape, fill = Cluster)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ Mutation, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = custom_colors_14) +
  labs(
    x = "Cluster",
    y = "Normalized Escape Score",
    title = "Escape scores by cluster for each mutation appearing in lineages of interest"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave("escape_boxplots_by_cluster.png", width = 12, height = 6, dpi = 300)
