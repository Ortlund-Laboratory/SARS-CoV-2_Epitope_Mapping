library(tidyverse)
library(dplyr)

dat <- read_csv("raw_data_cluster_app_over_time.csv")

# Wide matrix: rows = clusters, cols = months, values = % of cluster total
mat <- dat %>%
  mutate(
    Month   = as.integer(Month),
    Cluster = factor(Cluster)
  ) %>%
  select(Month, Cluster, Percentage_of_Cluster_Total) %>%
  pivot_wider(
    names_from = Month,
    values_from = Percentage_of_Cluster_Total
  ) %>%
  arrange(as.integer(Cluster)) %>%
  column_to_rownames("Cluster") %>%
  as.matrix()

# Build a label matrix with formatted percentages, e.g. "37%"
lab_mat <- apply(mat, 2, function(x) sprintf("%.0f", x))

library(pheatmap)

pheatmap(
  mat,
  cluster_rows = TRUE,           # Dendrogram
  cluster_cols = FALSE,          # Fixed month order
  breaks = seq(0, 100, length.out = 101),
  color = colorRampPalette(c("#FFFDF5", "#DAB1DA"))(100),
  #main = "Cluster representation over time",
  angle_col = 90,
  # show numbers in each cell
  display_numbers = lab_mat,
  number_color = "black",
  fontsize_number = 10,   # adjust to taste
  fontsize_col = 18,
  filename = "cluster_time_heatmap.png",  # ← Saves directly to PNG
  width = 12,
  height = 8
)

#Some additional file manipulations
# I want to determine a metric for each cluster as to the average number of antibodies
# studied over the 37 months

dat_cluster1 <- dat %>%
  filter(Cluster == 1)
mean_cluster1 <- mean(dat_cluster1$Percentage_of_Cluster_Total, na.rm = TRUE)

dat_cluster2 <- dat %>%
  filter(Cluster == 2)
mean_cluster2 <- mean(dat_cluster2$Percentage_of_Cluster_Total, na.rm = TRUE)

dat_cluster3 <- dat %>%
  filter(Cluster == 3)
mean_cluster3 <- mean(dat_cluster3$Percentage_of_Cluster_Total, na.rm = TRUE)

dat_cluster4 <- dat %>%
  filter(Cluster == 4)
mean_cluster4 <- mean(dat_cluster4$Percentage_of_Cluster_Total, na.rm = TRUE)

dat_cluster5 <- dat %>%
  filter(Cluster == 5)
mean_cluster5 <- mean(dat_cluster5$Percentage_of_Cluster_Total, na.rm = TRUE)

dat_cluster6 <- dat %>%
  filter(Cluster == 6)
mean_cluster6 <- mean(dat_cluster6$Percentage_of_Cluster_Total, na.rm = TRUE)

dat_cluster7 <- dat %>%
  filter(Cluster == 7)
mean_cluster7 <- mean(dat_cluster7$Percentage_of_Cluster_Total, na.rm = TRUE)

dat_cluster8 <- dat %>%
  filter(Cluster == 8)
mean_cluster8 <- mean(dat_cluster8$Percentage_of_Cluster_Total, na.rm = TRUE)

dat_cluster9 <- dat %>%
  filter(Cluster == 9)
mean_cluster9 <- mean(dat_cluster9$Percentage_of_Cluster_Total, na.rm = TRUE)

dat_cluster10 <- dat %>%
  filter(Cluster == 10)
mean_cluster10 <- mean(dat_cluster10$Percentage_of_Cluster_Total, na.rm = TRUE)

dat_cluster11 <- dat %>%
  filter(Cluster == 11)
mean_cluster11 <- mean(dat_cluster11$Percentage_of_Cluster_Total, na.rm = TRUE)

dat_cluster12 <- dat %>%
  filter(Cluster == 12)
mean_cluster12 <- mean(dat_cluster12$Percentage_of_Cluster_Total, na.rm = TRUE)

dat_cluster13 <- dat %>%
  filter(Cluster == 13)
mean_cluster13 <- mean(dat_cluster13$Percentage_of_Cluster_Total, na.rm = TRUE)

dat_cluster14 <- dat %>%
  filter(Cluster == 14)
mean_cluster14 <- mean(dat_cluster14$Percentage_of_Cluster_Total, na.rm = TRUE)
