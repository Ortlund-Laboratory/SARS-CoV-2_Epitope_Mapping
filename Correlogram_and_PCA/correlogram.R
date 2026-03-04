#install.packages("corrgram")
library(corrgram)

# set to your folder
path <- "/home/adkeith@Eu.Emory.Edu/DMS_Workflow/nucleocapsid/Covid_Update_Correlogram/cleanup_correlogram"

# list only files ending with *_escape_fractions_raw.csv
files <- list.files(
  path = path,
  pattern = "_escape_fractions_raw\\.csv$",
  full.names = TRUE
)

# derive names from the * part
names_vec <- sub("_escape_fractions_raw\\.csv$", "", basename(files))

# read all files into a named list
dat_list <- lapply(files, function(f) {
  read.csv(f, stringsAsFactors = FALSE)
})

# name each element after the * part
names(dat_list) <- names_vec

# OPTIONAL: if you really want separate objects in the global env
# (less recommended than keeping a list)
for (i in seq_along(dat_list)) {
  assign(names_vec[i], dat_list[[i]])
}

#Read in reference which contains all possible mutations
reference <- read.csv("reference.csv", header=TRUE)

##### Now combine dataframes into one
# put all your mutation data frames in a list
df_list <- mget(names_vec)   # etc.

# df_list is a *named* list; names(df_list) should match names_vec
# e.g. names(df_list) <- names_vec (if not already)
names(df_list) <- names_vec

for (nm in names(df_list)) {
  # rename Escape_Score -> <datasetname>_Escape_Score
  names(df_list[[nm]])[names(df_list[[nm]]) == "Escape_Score"] <-
    paste0(nm)
}

# Add row ID to preserve order
reference$row_id <- 1:nrow(reference)
combined <- reference

for (d in df_list) {
  combined <- merge(
    combined, 
    d, 
    by = "Mutation", 
    all.x = TRUE, 
    all.y = FALSE,
    sort = FALSE
  )
}

# Restore EXACT reference order
combined <- combined[order(combined$row_id), ]
combined$row_id <- NULL
rownames(combined) <- NULL

# Keep only numeric columns (corrgram requires this)
combined_num <- combined[, sapply(combined, is.numeric)]

### NEED TO REMOVE SITES 51-53, 71-73, 84-87, 110-116, 130 and 132
### In combined_num, this corresponds to lines 951-1007, 1331-1387, 1578-1653, 
### 2072-2204, 2452-2470 and 2490-2508.

remove_idx <- c(951:1007, 1331:1387, 1578:1653, 2072:2204, 2452:2470, 2490:2508)
combined_num_rem <- combined_num[-remove_idx, ]

# reorder
order_df   <- read.csv("corr_order.csv", header = TRUE)
order_vec  <- order_df$name      # character vector of variable names

combined_num_ord <- combined_num_rem[ , order_vec]

png("corrgram_final.png", width = 3000, height = 3000, res = 200)

par(oma = c(6, 6, 2, 6),    # more room right
    mar = c(2, 2, 4, 2),
    xpd = NA)

corrgram(
  combined_num_ord,
  order       = FALSE,
  lower.panel = panel.shade,
  upper.panel = NULL,
  labels      = colnames(combined_num_ord),
  label.pos   = c(2.4, 3.6),   # keep inside cell, near upper‑right
  cex.labels  = 1,
  label.srt   = 45,
  main        = "Antibody Correlation Matrix"
)

dev.off()



## Now create matrix with Pearson numbers

# Pearson correlation matrix (same values used by corrgram)
cmat <- cor(combined_num_ord, use = "pairwise.complete.obs", method = "pearson")

# round to 3 decimal places
cmat_round <- round(cmat, 3)

# look at it
cmat_round[1:5, 1:5]

write.csv(cmat_round, file = "pearson_correlations_round3.csv")

# cmat is your Pearson correlation matrix from: cmat <- cor(combined_num, use="pairwise.complete.obs")

# put correlations in a data frame, remove self‑correlations and duplicates
df_cor <- as.data.frame(as.table(cmat))
names(df_cor) <- c("var1", "var2", "r")

# drop diagonal and keep only one triangle
df_cor <- df_cor[df_cor$var1 != df_cor$var2, ]
df_cor <- df_cor[as.character(df_cor$var1) < as.character(df_cor$var2), ]

# order by absolute correlation, descending
df_cor <- df_cor[order(-abs(df_cor$r)), ]

# Top 50 by r
top50 <- df_cor[order(-df_cor$r), ][1:50, ]
top50$r <- round(top50$r, 3)

# Bottom 50 by r
bottom50 <- df_cor[order(df_cor$r), ][1:50, ]
bottom50$r <- round(bottom50$r, 3)

top50
bottom50


### Now let's look at PER-SITE correlations

n <- nrow(combined_num_ord)
block_size <- 19

# group index: 1 for rows 1–19, 2 for 20–38, etc.
g <- ((seq_len(n) - 1) %/% block_size) + 1

# compute column means within each block
avg_mat <- sapply(split(seq_len(n), g), function(idx) {
  colMeans(combined_num_ord[idx, , drop = FALSE], na.rm = TRUE)
})

# sapply gives variables in rows and groups in columns; transpose to keep
# same orientation as original: one row per 19-row block
avg_mat <- t(avg_mat)

# Remove rows that are entirely NA
avg_mat_clean <- avg_mat[!rowSums(is.na(avg_mat)) == ncol(avg_mat), ]

png("corrgram_final_per_site.png", width = 3000, height = 3000, res = 200)

par(oma = c(6, 6, 2, 6),    # more room right
    mar = c(2, 2, 4, 2),
    xpd = NA)

corrgram(
  avg_mat_clean,
  order       = FALSE,
  lower.panel = panel.shade,
  upper.panel = NULL,
  labels      = colnames(avg_mat_clean),
  label.pos   = c(2.4, 3.6),   # keep inside cell, near upper‑right
  cex.labels  = 1,
  label.srt   = 45,
  main        = "Per-Site Antibody Correlation Matrix"
)

dev.off()



## Now create matrix with Pearson numbers

# Pearson correlation matrix (same values used by corrgram)
cmatps <- cor(avg_mat_clean, use = "pairwise.complete.obs", method = "pearson")

# round to 3 decimal places
cmat_roundps <- round(cmatps, 3)

# look at it
cmat_roundps[1:5, 1:5]

write.csv(cmat_roundps, file = "per_site_pearson_correlations_round3.csv")

# put correlations in a data frame, remove self‑correlations and duplicates
df_corps <- as.data.frame(as.table(cmatps))
names(df_corps) <- c("var1", "var2", "r")

# drop diagonal and keep only one triangle
df_corps <- df_corps[df_corps$var1 != df_corps$var2, ]
df_corps <- df_corps[as.character(df_corps$var1) < as.character(df_corps$var2), ]

# order by absolute correlation, descending
df_corps <- df_corps[order(-abs(df_corps$r)), ]

# top 50, rounded to 3 decimals
top50ps <- head(df_corps, 50)
top50ps$r <- round(top50ps$r, 3)

top50ps


#### PCA for PER-VARIANT dataset

#### IT'S VERY IMPORTANT THAT ALL DATASETS ARE Z-NORMALIZED PRIOR TO PCA
#### OTHERWISE ANALYSIS WILL BE SEVERELY BIASED

# Z-normalize ALL columns (sites across all antibodies)
combined_num_ord_normalized <- as.data.frame(scale(combined_num_ord))

# Check it worked (means ~0, SDs ~1 per column)
colMeans(combined_num_ord_normalized, na.rm=TRUE)  # Should be ~0
apply(combined_num_ord_normalized, 2, sd, na.rm=TRUE)  # Should be ~1


# Transpose so antibodies are rows, sites are columns
data_pca <- t(combined_num_ord_normalized)  # sites x antibodies becomes antibodies x sites

# 1. Remove constant columns
keep_cols <- apply(data_pca, 2, function(x) var(x, na.rm=TRUE) > 0)
data_pca <- data_pca[, keep_cols]

# 2. Replace NAs with 0 (or median)
data_pca[is.na(data_pca)] <- 0

# 3. Remove infinite values (rows and columns)
finite_rows <- !rowSums(is.infinite(data_pca))
finite_cols <- !colSums(is.infinite(data_pca))
data_clean <- data_pca[finite_rows, finite_cols]

# 4. Final PCA
pca_res <- prcomp(data_clean, center = TRUE, scale. = FALSE)

summary(pca_res)

library(ggplot2)

# PC1 vs PC2 colored by cluster
plot_data <- data.frame(
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  antibody = rownames(data_pca)
)

# Initial clustering to see groups
k <- 5  # adjust based on elbow plot
clusters <- kmeans(pca_res$x[,1:3], centers = k)$cluster

ggplot(plot_data, aes(PC1, PC2, color = factor(clusters))) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Antibody epitope classes by PCA")

# Create metadata from PCA rownames (your antibody names)
antibody_metadata <- data.frame(
  antibody = rownames(pca_res$x),  # antibody names from PCA
  cluster = NA  # will fill this
)

#Optimal clusters, elbow method
library(cluster)

# Use first 4 PCs (explains most epitope variance)
pca_data <- pca_res$x[,1:4]

# Elbow method
wss <- sapply(1:50, function(k) {
  kmeans(pca_data, centers=k, nstart=25)$tot.withinss
})

# Plot elbow
plot_data <- data.frame(k=1:50, wss=wss)
ggplot(plot_data, aes(k, wss)) + 
  geom_point(size=3) + geom_line() +
  labs(x="Number of Clusters", y="Within-cluster SS") +
  theme_minimal(base_size=14)
# Look for elbow (usually 3-5 for epitopes)
# Yes, this looks correct. But let's see what sillhouette suggests

# Normalized WSS (1 = k=1, approaching 0 = many clusters)
wss_norm <- wss / max(wss)

# Find where 80% of "elbow drop" has occurred
elbow_80 <- min(which(wss_norm <= 0.2))
cat("80% elbow rule: k =", elbow_80, "\n")

ggsave("elbow.png")

#With a fixed seed and reasonably large nstart / number of repeats, your 
#silhouette curve should become highly reproducible

set.seed(123)
sil <- sapply(2:60, function(k) {
  scores <- replicate(20, {
    km <- kmeans(pca_data, centers = k, nstart = 10)
    mean(silhouette(km$cluster, dist(pca_data))[, 3])
  })
  max(scores)   # or mean(scores)
})

# Save silhouette plot
png("silhouette_plot.png", width=8, height=6, units="in", res=300, pointsize=14)

plot(2:60, sil, type="b", pch=19, xlab="Number of Clusters (k)", 
     ylab="Average Silhouette Width", main="Optimal Cluster Number")
abline(v=which.max(sil)+1, col="red", lty=2, lwd=2)
text(which.max(sil)+1, max(sil), paste("Optimal k =", 22), pos=4, col="red")
#text(which.max(sil)+1, max(sil), paste("Optimal k ="), pos=4, col="red")

dev.off()
cat("Saved silhouette_plot.png\n")

# Silhouette suggests 22. But that's mainly noise
# Elbow method suggests 8

# Let's try gap statistic:
library(cluster)

# After running clusGap()
gap_stat <- clusGap(pca_data, FUN = kmeans, nstart = 25, K.max = 20, B = 50)

# Extract optimal k (automatic)
optimal_k <- maxSE(gap_stat$Tab[,"gap"], gap_stat$Tab[,"SE.sim"], method="firstmax")
cat("Optimal k =", optimal_k, "\n")

#dev.off()

# Gap method optimal k is 14.
# So box plot method gives 8
# Silhouette method gives 22
# Gap method is 14
# Gap method seems to have the most statistical rigor, plus it's roughly
# intermediate between the two other methods.
# So let's set optimal_k as 14

optimal_k <- 14

best_sil <- -Inf
best_km <- NULL

for (i in 1:20) {
  km <- kmeans(pca_data, centers = 14, nstart = 100)
  sil <- silhouette(km$cluster, dist(pca_data))
  avg_sil_i <- mean(sil[, 3])
  if (avg_sil_i > best_sil) {
    best_sil <- avg_sil_i
    best_km <- km
  }
}

pca_data_3pc <- pca_res$x[,1:3]
kmeans_3pc <- kmeans(pca_data_3pc, 14, nstart=100)
sil_3pc <- silhouette(kmeans_3pc$cluster, dist(pca_data_3pc))
mean(sil_3pc[,3])  # Might improve from 0.41
# Yes, it goes up to 0.46
# Your new gold standard
pca_data_3pc <- pca_res$x[,1:3]  # PC1+PC2+PC3
final_clusters <- kmeans(pca_data_3pc, centers=14, nstart=100)
antibody_metadata$cluster <- final_clusters$cluster

# Confirm the win
sil_final <- silhouette(final_clusters$cluster, dist(pca_data_3pc))
cat("Final avg silhouette:", round(mean(sil_final[,3]), 3), "↑ from 0.406\n")

pca_df_3d <- data.frame(PC1=pca_res$x[,1], PC2=pca_res$x[,2], PC3=pca_res$x[,3],
                        Cluster=factor(final_clusters$cluster))

# Custom color palette (cluster 6 = dark teal instead of yellow)
custom_colors_14 <- c("#4DAF4A", "#E4C899", "#377EB8", "#E41A1C", 
                      "#A65628", "#984EA3", "grey40", "#D3AF37",
                      # Extensions (complementary tones)
                      "#F781BF", "#999999", "#66C2A5", 
                      "#FF7F00", "#b1f2ff", "#03045E")

ggplot(pca_df_3d, aes(PC1, PC2, color=Cluster)) + 
  geom_point(size=2.5, alpha=0.75) +
  scale_color_manual(values=custom_colors_14) +
  labs(title="Epitope Clusters (k=14, 3-PC optimized, sil=0.46)",
       subtitle="PC1:14% + PC2:5% + PC3 improved separation") +
  theme_minimal(base_size=14)

library(ggrepel)

# Create 3D PCA data frame with labels
pca_df_3d <- data.frame(
  PC1 = pca_res$x[,1], 
  PC2 = pca_res$x[,2], 
  PC3 = pca_res$x[,3],
  Cluster = factor(final_clusters$cluster),
  Antibody = antibody_metadata$antibody
)

# 2D plot of PC1 vs PC2 with 3D clustering (labels all antibodies)
p_labeled <- ggplot(pca_df_3d, aes(PC1, PC2, color=Cluster)) + 
  geom_point(size=3, alpha=0.8) +
  geom_text_repel(aes(label=Antibody), 
                  size=3, 
                  max.overlaps=20,
                  segment.size=0.2,
                  segment.alpha=0.5,
                  box.padding=0.3) +
  scale_color_manual(values=custom_colors_14) +
  labs(x=paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100,1), "%)"),
       y=paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100,1), "%)"),
       title="Epitope Clusters (k=14, 3-PC optimized, sil=0.46)",
       subtitle="All antibodies labeled") +
  theme_minimal(base_size=14) +
  theme(legend.position="bottom")

# Save high-res PNG
ggsave("epitope_pca_3pc_labeled.png", 
       plot=p_labeled, 
       width=16, height=12, 
       dpi=300, 
       bg="white")


########################################################
#THIS COULD BE USEFUL BUT WILL NEED TO GIVE SOME THOUGHT
########################################################

#Find top escape sites per cluster

## AFTER avg_mat_clean is created
## Build antibody × site matrix

# I want to use z-normalized data
n <- nrow(combined_num_ord_normalized)
block_size <- 19

# group index: 1 for rows 1–19, 2 for 20–38, etc.
g <- ((seq_len(n) - 1) %/% block_size) + 1

# compute column means within each block
avg_mat_norm <- sapply(split(seq_len(n), g), function(idx) {
  colMeans(combined_num_ord_normalized[idx, , drop = FALSE], na.rm = TRUE)
})

# sapply gives variables in rows and groups in columns; transpose to keep
# same orientation as original: one row per 19-row block
avg_mat_norm <- t(avg_mat_norm)

# Remove rows that are entirely NA
avg_mat_clean_norm <- avg_mat_norm[!rowSums(is.na(avg_mat_norm)) == ncol(avg_mat_norm), ]

# avg_mat_clean: rows = sites, columns = antibodies
site_by_antibody <- t(avg_mat_clean_norm)   # rows = antibodies, cols = sites

# Get current numbers
nums <- as.numeric(sub("Site_", "", colnames(site_by_antibody)))

# Apply each gap correction sequentially
nums[nums >= 51] <- nums[nums >= 51] + 3   # Skip 51,52,53  
nums[nums >= 71] <- nums[nums >= 71] + 3   # Skip 71,72,73
nums[nums >= 84] <- nums[nums >= 84] + 4   # Skip 84-87  
nums[nums >= 110] <- nums[nums >= 110] + 7 # Skip 110-116
nums[nums >= 130] <- nums[nums >= 130] + 1 # Skip 130
nums[nums >= 132] <- nums[nums >= 132] + 1 # Skip 130

# Put back
colnames(site_by_antibody) <- paste0("Site_", nums)

# Check
head(colnames(site_by_antibody), 60)
tail(colnames(site_by_antibody))
max(nums)  # Should be 419



# Align antibodies with antibody_metadata
site_by_antibody <- site_by_antibody[match(antibody_metadata$antibody,
                                           rownames(site_by_antibody)), ]

stopifnot(all(rownames(site_by_antibody) == antibody_metadata$antibody))


escape_site <- site_by_antibody
clusters    <- antibody_metadata$cluster
site_names  <- colnames(escape_site)
site_names
top_n       <- 10

cluster_site_means_list <- lapply(sort(unique(clusters)), function(k) {
  idx   <- clusters == k
  mat_k <- escape_site[idx, , drop = FALSE]
  mean_k <- colMeans(mat_k, na.rm = TRUE)
  ord <- order(mean_k, decreasing = TRUE)[1:top_n]
  
  # FIXED: Create data frame properly with separate columns
  data.frame(
    Cluster     = rep(k, top_n),
    Site        = site_names[ord],
    Mean_escape = mean_k[ord]
  )
})

cluster_site_means_df <- do.call(rbind, cluster_site_means_list)
rownames(cluster_site_means_df) <- NULL  # Clean row names
cluster_site_means_df

# Two standard deviation above the mean
# for all the sites in each cluster, with a condition that these identified
# sites should also be at least 0.5SD above the mean for EACH member of the cluster.

# 1. Put escape_site into a data frame and attach clusters
escape_df <- as.data.frame(escape_site)
escape_df$Cluster <- clusters  # same order as rows

# Optional: keep company IDs as a column instead of rownames
escape_df$Company <- rownames(escape_site)

# 2. Split into a list of data frames by cluster
escape_by_cluster <- split(escape_df, escape_df$Cluster)

# Now I have:
# escape_by_cluster[["1"]]  -> all companies in cluster 1
# escape_by_cluster[["2"]]  -> all companies in cluster 2
# ...
# escape_by_cluster[["14"]] -> all companies in cluster 14

# Get just the escape matrix (drop Cluster/Company columns)
escape_mat_cluster_1 <- as.matrix(
  escape_by_cluster[["1"]][, !(names(escape_by_cluster[["1"]]) %in% c("Cluster","Company"))]
)
escape_mat_cluster_2 <- as.matrix(
  escape_by_cluster[["2"]][, !(names(escape_by_cluster[["2"]]) %in% c("Cluster","Company"))]
)
escape_mat_cluster_3 <- as.matrix(
  escape_by_cluster[["3"]][, !(names(escape_by_cluster[["3"]]) %in% c("Cluster","Company"))]
)
escape_mat_cluster_4 <- as.matrix(
  escape_by_cluster[["4"]][, !(names(escape_by_cluster[["4"]]) %in% c("Cluster","Company"))]
)
escape_mat_cluster_5 <- as.matrix(
  escape_by_cluster[["5"]][, !(names(escape_by_cluster[["5"]]) %in% c("Cluster","Company"))]
)
escape_mat_cluster_6 <- as.matrix(
  escape_by_cluster[["6"]][, !(names(escape_by_cluster[["6"]]) %in% c("Cluster","Company"))]
)
escape_mat_cluster_7 <- as.matrix(
  escape_by_cluster[["7"]][, !(names(escape_by_cluster[["7"]]) %in% c("Cluster","Company"))]
)
escape_mat_cluster_8 <- as.matrix(
  escape_by_cluster[["8"]][, !(names(escape_by_cluster[["8"]]) %in% c("Cluster","Company"))]
)
escape_mat_cluster_9 <- as.matrix(
  escape_by_cluster[["9"]][, !(names(escape_by_cluster[["9"]]) %in% c("Cluster","Company"))]
)
escape_mat_cluster_10 <- as.matrix(
  escape_by_cluster[["10"]][, !(names(escape_by_cluster[["10"]]) %in% c("Cluster","Company"))]
)
escape_mat_cluster_11 <- as.matrix(
  escape_by_cluster[["11"]][, !(names(escape_by_cluster[["11"]]) %in% c("Cluster","Company"))]
)
escape_mat_cluster_12 <- as.matrix(
  escape_by_cluster[["12"]][, !(names(escape_by_cluster[["12"]]) %in% c("Cluster","Company"))]
)
escape_mat_cluster_13 <- as.matrix(
  escape_by_cluster[["13"]][, !(names(escape_by_cluster[["13"]]) %in% c("Cluster","Company"))]
)
escape_mat_cluster_14 <- as.matrix(
  escape_by_cluster[["14"]][, !(names(escape_by_cluster[["14"]]) %in% c("Cluster","Company"))]
)



# Using escape_mat_cluster (matrix: rows = companies, cols = sites)
############## CLUSTER 1 ########################


# 1. Cluster-level: 2 SD threshold
mean_escape_per_site <- colMeans(escape_mat_cluster_1, na.rm = TRUE)
mean_escape_per_site <- mean_escape_per_site[!is.nan(mean_escape_per_site)]

overall_mean <- mean(mean_escape_per_site)
site_sd <- sd(mean_escape_per_site)
cluster_threshold <- overall_mean + 2 * site_sd

# LIST 1: Cluster-only (2SD threshold)
cluster_only_sites <- names(mean_escape_per_site)[mean_escape_per_site >= cluster_threshold]

# 2. Per-company stats (0.5 SD threshold)
per_company_stats <- lapply(1:nrow(escape_mat_cluster_1), function(i) {
  company_row <- escape_mat_cluster_1[i, , drop = FALSE]
  company_mean <- mean(company_row, na.rm = TRUE)
  company_sd <- sd(company_row, na.rm = TRUE)
  list(mean = company_mean, sd = company_sd)
})

# LIST 2: Universal sites (2SD cluster + 0.5SD EVERY company)
universal_sig_sites <- cluster_only_sites
for(i in 1:nrow(escape_mat_cluster_1)) {
  company_row <- as.numeric(escape_mat_cluster_1[i, cluster_only_sites])
  company_threshold <- per_company_stats[[i]]$mean + 0.5 * per_company_stats[[i]]$sd
  
  universal_sig_sites <- universal_sig_sites[company_row >= company_threshold]
}

# Summary
cat("=== CLUSTER 1 EPITOPE SITES ===\n")
cat("Cluster threshold (2 SD):", round(cluster_threshold, 4), "\n")
cat("\nLIST 1 - Cluster-only (2SD):", length(cluster_only_sites), "sites\n")
print(sort(cluster_only_sites))

cat("\nLIST 2 - Universal (2SD cluster + 0.5SD all companies):", 
    length(universal_sig_sites), "sites\n")
print(sort(universal_sig_sites))




# Using escape_mat_cluster (matrix: rows = companies, cols = sites)
############## CLUSTER 2 ########################


# 1. Cluster-level: 2 SD threshold
mean_escape_per_site <- colMeans(escape_mat_cluster_2, na.rm = TRUE)
mean_escape_per_site <- mean_escape_per_site[!is.nan(mean_escape_per_site)]

overall_mean <- mean(mean_escape_per_site)
site_sd <- sd(mean_escape_per_site)
cluster_threshold <- overall_mean + 2 * site_sd

# LIST 1: Cluster-only (2SD threshold)
cluster_only_sites <- names(mean_escape_per_site)[mean_escape_per_site >= cluster_threshold]

# 2. Per-company stats (0.5 SD threshold)
per_company_stats <- lapply(1:nrow(escape_mat_cluster_2), function(i) {
  company_row <- escape_mat_cluster_2[i, , drop = FALSE]
  company_mean <- mean(company_row, na.rm = TRUE)
  company_sd <- sd(company_row, na.rm = TRUE)
  list(mean = company_mean, sd = company_sd)
})

# LIST 2: Universal sites (2SD cluster + 0.5SD EVERY company)
universal_sig_sites <- cluster_only_sites
for(i in 1:nrow(escape_mat_cluster_2)) {
  company_row <- as.numeric(escape_mat_cluster_2[i, cluster_only_sites])
  company_threshold <- per_company_stats[[i]]$mean + 0.5 * per_company_stats[[i]]$sd
  
  universal_sig_sites <- universal_sig_sites[company_row >= company_threshold]
}

# Summary
cat("=== CLUSTER 2 EPITOPE SITES ===\n")
cat("Cluster threshold (2 SD):", round(cluster_threshold, 4), "\n")
cat("\nLIST 1 - Cluster-only (2SD):", length(cluster_only_sites), "sites\n")
print(sort(cluster_only_sites))

cat("\nLIST 2 - Universal (2SD cluster + 0.5SD all companies):", 
    length(universal_sig_sites), "sites\n")
print(sort(universal_sig_sites))

# Using escape_mat_cluster (matrix: rows = companies, cols = sites)
############## CLUSTER 3 ########################


# 1. Cluster-level: 2 SD threshold
mean_escape_per_site <- colMeans(escape_mat_cluster_3, na.rm = TRUE)
mean_escape_per_site <- mean_escape_per_site[!is.nan(mean_escape_per_site)]

overall_mean <- mean(mean_escape_per_site)
site_sd <- sd(mean_escape_per_site)
cluster_threshold <- overall_mean + 2 * site_sd

# LIST 1: Cluster-only (2SD threshold)
cluster_only_sites <- names(mean_escape_per_site)[mean_escape_per_site >= cluster_threshold]

# 2. Per-company stats (0.5 SD threshold)
per_company_stats <- lapply(1:nrow(escape_mat_cluster_3), function(i) {
  company_row <- escape_mat_cluster_3[i, , drop = FALSE]
  company_mean <- mean(company_row, na.rm = TRUE)
  company_sd <- sd(company_row, na.rm = TRUE)
  list(mean = company_mean, sd = company_sd)
})

# LIST 2: Universal sites (2SD cluster + 0.5SD EVERY company)
universal_sig_sites <- cluster_only_sites
for(i in 1:nrow(escape_mat_cluster_3)) {
  company_row <- as.numeric(escape_mat_cluster_3[i, cluster_only_sites])
  company_threshold <- per_company_stats[[i]]$mean + 0.5 * per_company_stats[[i]]$sd
  
  universal_sig_sites <- universal_sig_sites[company_row >= company_threshold]
}

# Summary
cat("=== CLUSTER 3 EPITOPE SITES ===\n")
cat("Cluster threshold (2 SD):", round(cluster_threshold, 4), "\n")
cat("\nLIST 1 - Cluster-only (2SD):", length(cluster_only_sites), "sites\n")
print(sort(cluster_only_sites))

cat("\nLIST 2 - Universal (2SD cluster + 0.5SD all companies):", 
    length(universal_sig_sites), "sites\n")
print(sort(universal_sig_sites))

# Using escape_mat_cluster (matrix: rows = companies, cols = sites)
############## CLUSTER 4 ########################


# 1. Cluster-level: 2 SD threshold
mean_escape_per_site <- colMeans(escape_mat_cluster_4, na.rm = TRUE)
mean_escape_per_site <- mean_escape_per_site[!is.nan(mean_escape_per_site)]

overall_mean <- mean(mean_escape_per_site)
site_sd <- sd(mean_escape_per_site)
cluster_threshold <- overall_mean + 2 * site_sd

# LIST 1: Cluster-only (2SD threshold)
cluster_only_sites <- names(mean_escape_per_site)[mean_escape_per_site >= cluster_threshold]

# 2. Per-company stats (0.5 SD threshold)
per_company_stats <- lapply(1:nrow(escape_mat_cluster_4), function(i) {
  company_row <- escape_mat_cluster_4[i, , drop = FALSE]
  company_mean <- mean(company_row, na.rm = TRUE)
  company_sd <- sd(company_row, na.rm = TRUE)
  list(mean = company_mean, sd = company_sd)
})

# LIST 2: Universal sites (2SD cluster + 0.5SD EVERY company)
universal_sig_sites <- cluster_only_sites
for(i in 1:nrow(escape_mat_cluster_4)) {
  company_row <- as.numeric(escape_mat_cluster_4[i, cluster_only_sites])
  company_threshold <- per_company_stats[[i]]$mean + 0.5 * per_company_stats[[i]]$sd
  
  universal_sig_sites <- universal_sig_sites[company_row >= company_threshold]
}

# Summary
cat("=== CLUSTER 4 EPITOPE SITES ===\n")
cat("Cluster threshold (2 SD):", round(cluster_threshold, 4), "\n")
cat("\nLIST 1 - Cluster-only (2SD):", length(cluster_only_sites), "sites\n")
print(sort(cluster_only_sites))

cat("\nLIST 2 - Universal (2SD cluster + 0.5SD all companies):", 
    length(universal_sig_sites), "sites\n")
print(sort(universal_sig_sites))

# Using escape_mat_cluster (matrix: rows = companies, cols = sites)
############## CLUSTER 5 ########################


# 1. Cluster-level: 2 SD threshold
mean_escape_per_site <- colMeans(escape_mat_cluster_5, na.rm = TRUE)
mean_escape_per_site <- mean_escape_per_site[!is.nan(mean_escape_per_site)]

overall_mean <- mean(mean_escape_per_site)
site_sd <- sd(mean_escape_per_site)
cluster_threshold <- overall_mean + 2 * site_sd

# LIST 1: Cluster-only (2SD threshold)
cluster_only_sites <- names(mean_escape_per_site)[mean_escape_per_site >= cluster_threshold]

# 2. Per-company stats (0.5 SD threshold)
per_company_stats <- lapply(1:nrow(escape_mat_cluster_5), function(i) {
  company_row <- escape_mat_cluster_5[i, , drop = FALSE]
  company_mean <- mean(company_row, na.rm = TRUE)
  company_sd <- sd(company_row, na.rm = TRUE)
  list(mean = company_mean, sd = company_sd)
})

# LIST 2: Universal sites (2SD cluster + 0.5SD EVERY company)
universal_sig_sites <- cluster_only_sites
for(i in 1:nrow(escape_mat_cluster_5)) {
  company_row <- as.numeric(escape_mat_cluster_5[i, cluster_only_sites])
  company_threshold <- per_company_stats[[i]]$mean + 0.5 * per_company_stats[[i]]$sd
  
  universal_sig_sites <- universal_sig_sites[company_row >= company_threshold]
}

# Summary
cat("=== CLUSTER 5 EPITOPE SITES ===\n")
cat("Cluster threshold (2 SD):", round(cluster_threshold, 4), "\n")
cat("\nLIST 1 - Cluster-only (2SD):", length(cluster_only_sites), "sites\n")
print(sort(cluster_only_sites))

cat("\nLIST 2 - Universal (2SD cluster + 0.5SD all companies):", 
    length(universal_sig_sites), "sites\n")
print(sort(universal_sig_sites))

# Using escape_mat_cluster (matrix: rows = companies, cols = sites)
############## CLUSTER 6 ########################


# 1. Cluster-level: 2 SD threshold
mean_escape_per_site <- colMeans(escape_mat_cluster_6, na.rm = TRUE)
mean_escape_per_site <- mean_escape_per_site[!is.nan(mean_escape_per_site)]

overall_mean <- mean(mean_escape_per_site)
site_sd <- sd(mean_escape_per_site)
cluster_threshold <- overall_mean + 2 * site_sd

# LIST 1: Cluster-only (2SD threshold)
cluster_only_sites <- names(mean_escape_per_site)[mean_escape_per_site >= cluster_threshold]

# 2. Per-company stats (0.5 SD threshold)
per_company_stats <- lapply(1:nrow(escape_mat_cluster_6), function(i) {
  company_row <- escape_mat_cluster_6[i, , drop = FALSE]
  company_mean <- mean(company_row, na.rm = TRUE)
  company_sd <- sd(company_row, na.rm = TRUE)
  list(mean = company_mean, sd = company_sd)
})

# LIST 2: Universal sites (2SD cluster + 0.5SD EVERY company)
universal_sig_sites <- cluster_only_sites
for(i in 1:nrow(escape_mat_cluster_6)) {
  company_row <- as.numeric(escape_mat_cluster_6[i, cluster_only_sites])
  company_threshold <- per_company_stats[[i]]$mean + 0.5 * per_company_stats[[i]]$sd
  
  universal_sig_sites <- universal_sig_sites[company_row >= company_threshold]
}

# Summary
cat("=== CLUSTER 6 EPITOPE SITES ===\n")
cat("Cluster threshold (2 SD):", round(cluster_threshold, 4), "\n")
cat("\nLIST 1 - Cluster-only (2SD):", length(cluster_only_sites), "sites\n")
print(sort(cluster_only_sites))

cat("\nLIST 2 - Universal (2SD cluster + 0.5SD all companies):", 
    length(universal_sig_sites), "sites\n")
print(sort(universal_sig_sites))

# Using escape_mat_cluster (matrix: rows = companies, cols = sites)
############## CLUSTER 7 ########################


# 1. Cluster-level: 2 SD threshold
mean_escape_per_site <- colMeans(escape_mat_cluster_7, na.rm = TRUE)
mean_escape_per_site <- mean_escape_per_site[!is.nan(mean_escape_per_site)]

overall_mean <- mean(mean_escape_per_site)
site_sd <- sd(mean_escape_per_site)
cluster_threshold <- overall_mean + 2 * site_sd

# LIST 1: Cluster-only (2SD threshold)
cluster_only_sites <- names(mean_escape_per_site)[mean_escape_per_site >= cluster_threshold]

# 2. Per-company stats (0.5 SD threshold)
per_company_stats <- lapply(1:nrow(escape_mat_cluster_7), function(i) {
  company_row <- escape_mat_cluster_7[i, , drop = FALSE]
  company_mean <- mean(company_row, na.rm = TRUE)
  company_sd <- sd(company_row, na.rm = TRUE)
  list(mean = company_mean, sd = company_sd)
})

# LIST 2: Universal sites (2SD cluster + 0.5SD EVERY company)
universal_sig_sites <- cluster_only_sites
for(i in 1:nrow(escape_mat_cluster_7)) {
  company_row <- as.numeric(escape_mat_cluster_7[i, cluster_only_sites])
  company_threshold <- per_company_stats[[i]]$mean + 0.5 * per_company_stats[[i]]$sd
  
  universal_sig_sites <- universal_sig_sites[company_row >= company_threshold]
}

# Summary
cat("=== CLUSTER 7 EPITOPE SITES ===\n")
cat("Cluster threshold (2 SD):", round(cluster_threshold, 4), "\n")
cat("\nLIST 1 - Cluster-only (2SD):", length(cluster_only_sites), "sites\n")
print(sort(cluster_only_sites))

cat("\nLIST 2 - Universal (2SD cluster + 0.5SD all companies):", 
    length(universal_sig_sites), "sites\n")
print(sort(universal_sig_sites))

# Using escape_mat_cluster (matrix: rows = companies, cols = sites)
############## CLUSTER 8 ########################


# 1. Cluster-level: 2 SD threshold
mean_escape_per_site <- colMeans(escape_mat_cluster_8, na.rm = TRUE)
mean_escape_per_site <- mean_escape_per_site[!is.nan(mean_escape_per_site)]

overall_mean <- mean(mean_escape_per_site)
site_sd <- sd(mean_escape_per_site)
cluster_threshold <- overall_mean + 2 * site_sd

# LIST 1: Cluster-only (2SD threshold)
cluster_only_sites <- names(mean_escape_per_site)[mean_escape_per_site >= cluster_threshold]

# 2. Per-company stats (0.5 SD threshold)
per_company_stats <- lapply(1:nrow(escape_mat_cluster_8), function(i) {
  company_row <- escape_mat_cluster_8[i, , drop = FALSE]
  company_mean <- mean(company_row, na.rm = TRUE)
  company_sd <- sd(company_row, na.rm = TRUE)
  list(mean = company_mean, sd = company_sd)
})

# LIST 2: Universal sites (2SD cluster + 0.5SD EVERY company)
universal_sig_sites <- cluster_only_sites
for(i in 1:nrow(escape_mat_cluster_8)) {
  company_row <- as.numeric(escape_mat_cluster_8[i, cluster_only_sites])
  company_threshold <- per_company_stats[[i]]$mean + 0.5 * per_company_stats[[i]]$sd
  
  universal_sig_sites <- universal_sig_sites[company_row >= company_threshold]
}

# Summary
cat("=== CLUSTER 8 EPITOPE SITES ===\n")
cat("Cluster threshold (2 SD):", round(cluster_threshold, 4), "\n")
cat("\nLIST 1 - Cluster-only (2SD):", length(cluster_only_sites), "sites\n")
print(sort(cluster_only_sites))

cat("\nLIST 2 - Universal (2SD cluster + 0.5SD all companies):", 
    length(universal_sig_sites), "sites\n")
print(sort(universal_sig_sites))

# Using escape_mat_cluster (matrix: rows = companies, cols = sites)
############## CLUSTER 9 ########################


# 1. Cluster-level: 2 SD threshold
mean_escape_per_site <- colMeans(escape_mat_cluster_9, na.rm = TRUE)
mean_escape_per_site <- mean_escape_per_site[!is.nan(mean_escape_per_site)]

overall_mean <- mean(mean_escape_per_site)
site_sd <- sd(mean_escape_per_site)
cluster_threshold <- overall_mean + 2 * site_sd

# LIST 1: Cluster-only (2SD threshold)
cluster_only_sites <- names(mean_escape_per_site)[mean_escape_per_site >= cluster_threshold]

# 2. Per-company stats (0.5 SD threshold)
per_company_stats <- lapply(1:nrow(escape_mat_cluster_9), function(i) {
  company_row <- escape_mat_cluster_9[i, , drop = FALSE]
  company_mean <- mean(company_row, na.rm = TRUE)
  company_sd <- sd(company_row, na.rm = TRUE)
  list(mean = company_mean, sd = company_sd)
})

# LIST 2: Universal sites (2SD cluster + 0.5SD EVERY company)
universal_sig_sites <- cluster_only_sites
for(i in 1:nrow(escape_mat_cluster_9)) {
  company_row <- as.numeric(escape_mat_cluster_9[i, cluster_only_sites])
  company_threshold <- per_company_stats[[i]]$mean + 0.5 * per_company_stats[[i]]$sd
  
  universal_sig_sites <- universal_sig_sites[company_row >= company_threshold]
}

# Summary
cat("=== CLUSTER 9 EPITOPE SITES ===\n")
cat("Cluster threshold (2 SD):", round(cluster_threshold, 4), "\n")
cat("\nLIST 1 - Cluster-only (2SD):", length(cluster_only_sites), "sites\n")
print(sort(cluster_only_sites))

cat("\nLIST 2 - Universal (2SD cluster + 0.5SD all companies):", 
    length(universal_sig_sites), "sites\n")
print(sort(universal_sig_sites))

# Using escape_mat_cluster (matrix: rows = companies, cols = sites)
############## CLUSTER 10 ########################


# 1. Cluster-level: 2 SD threshold
mean_escape_per_site <- colMeans(escape_mat_cluster_10, na.rm = TRUE)
mean_escape_per_site <- mean_escape_per_site[!is.nan(mean_escape_per_site)]

overall_mean <- mean(mean_escape_per_site)
site_sd <- sd(mean_escape_per_site)
cluster_threshold <- overall_mean + 2 * site_sd

# LIST 1: Cluster-only (2SD threshold)
cluster_only_sites <- names(mean_escape_per_site)[mean_escape_per_site >= cluster_threshold]

# 2. Per-company stats (0.5 SD threshold)
per_company_stats <- lapply(1:nrow(escape_mat_cluster_10), function(i) {
  company_row <- escape_mat_cluster_10[i, , drop = FALSE]
  company_mean <- mean(company_row, na.rm = TRUE)
  company_sd <- sd(company_row, na.rm = TRUE)
  list(mean = company_mean, sd = company_sd)
})

# LIST 2: Universal sites (2SD cluster + 0.5SD EVERY company)
universal_sig_sites <- cluster_only_sites
for(i in 1:nrow(escape_mat_cluster_10)) {
  company_row <- as.numeric(escape_mat_cluster_10[i, cluster_only_sites])
  company_threshold <- per_company_stats[[i]]$mean + 0.5 * per_company_stats[[i]]$sd
  
  universal_sig_sites <- universal_sig_sites[company_row >= company_threshold]
}

# Summary
cat("=== CLUSTER 10 EPITOPE SITES ===\n")
cat("Cluster threshold (2 SD):", round(cluster_threshold, 4), "\n")
cat("\nLIST 1 - Cluster-only (2SD):", length(cluster_only_sites), "sites\n")
print(sort(cluster_only_sites))

cat("\nLIST 2 - Universal (2SD cluster + 0.5SD all companies):", 
    length(universal_sig_sites), "sites\n")
print(sort(universal_sig_sites))

# Using escape_mat_cluster (matrix: rows = companies, cols = sites)
############## CLUSTER 11 ########################


# 1. Cluster-level: 2 SD threshold
mean_escape_per_site <- colMeans(escape_mat_cluster_11, na.rm = TRUE)
mean_escape_per_site <- mean_escape_per_site[!is.nan(mean_escape_per_site)]

overall_mean <- mean(mean_escape_per_site)
site_sd <- sd(mean_escape_per_site)
cluster_threshold <- overall_mean + 2 * site_sd

# LIST 1: Cluster-only (2SD threshold)
cluster_only_sites <- names(mean_escape_per_site)[mean_escape_per_site >= cluster_threshold]

# 2. Per-company stats (0.5 SD threshold)
per_company_stats <- lapply(1:nrow(escape_mat_cluster_11), function(i) {
  company_row <- escape_mat_cluster_11[i, , drop = FALSE]
  company_mean <- mean(company_row, na.rm = TRUE)
  company_sd <- sd(company_row, na.rm = TRUE)
  list(mean = company_mean, sd = company_sd)
})

# LIST 2: Universal sites (2SD cluster + 0.5SD EVERY company)
universal_sig_sites <- cluster_only_sites
for(i in 1:nrow(escape_mat_cluster_11)) {
  company_row <- as.numeric(escape_mat_cluster_11[i, cluster_only_sites])
  company_threshold <- per_company_stats[[i]]$mean + 0.5 * per_company_stats[[i]]$sd
  
  universal_sig_sites <- universal_sig_sites[company_row >= company_threshold]
}

# Summary
cat("=== CLUSTER 11 EPITOPE SITES ===\n")
cat("Cluster threshold (2 SD):", round(cluster_threshold, 4), "\n")
cat("\nLIST 1 - Cluster-only (2SD):", length(cluster_only_sites), "sites\n")
print(sort(cluster_only_sites))

cat("\nLIST 2 - Universal (2SD cluster + 0.5SD all companies):", 
    length(universal_sig_sites), "sites\n")
print(sort(universal_sig_sites))

# Using escape_mat_cluster (matrix: rows = companies, cols = sites)
############## CLUSTER 12 ########################


# 1. Cluster-level: 2 SD threshold
mean_escape_per_site <- colMeans(escape_mat_cluster_12, na.rm = TRUE)
mean_escape_per_site <- mean_escape_per_site[!is.nan(mean_escape_per_site)]

overall_mean <- mean(mean_escape_per_site)
site_sd <- sd(mean_escape_per_site)
cluster_threshold <- overall_mean + 2 * site_sd

# LIST 1: Cluster-only (2SD threshold)
cluster_only_sites <- names(mean_escape_per_site)[mean_escape_per_site >= cluster_threshold]

# 2. Per-company stats (0.5 SD threshold)
per_company_stats <- lapply(1:nrow(escape_mat_cluster_12), function(i) {
  company_row <- escape_mat_cluster_12[i, , drop = FALSE]
  company_mean <- mean(company_row, na.rm = TRUE)
  company_sd <- sd(company_row, na.rm = TRUE)
  list(mean = company_mean, sd = company_sd)
})

# LIST 2: Universal sites (2SD cluster + 0.5SD EVERY company)
universal_sig_sites <- cluster_only_sites
for(i in 1:nrow(escape_mat_cluster_12)) {
  company_row <- as.numeric(escape_mat_cluster_12[i, cluster_only_sites])
  company_threshold <- per_company_stats[[i]]$mean + 0.5 * per_company_stats[[i]]$sd
  
  universal_sig_sites <- universal_sig_sites[company_row >= company_threshold]
}

# Summary
cat("=== CLUSTER 12 EPITOPE SITES ===\n")
cat("Cluster threshold (2 SD):", round(cluster_threshold, 4), "\n")
cat("\nLIST 1 - Cluster-only (2SD):", length(cluster_only_sites), "sites\n")
print(sort(cluster_only_sites))

cat("\nLIST 2 - Universal (2SD cluster + 0.5SD all companies):", 
    length(universal_sig_sites), "sites\n")
print(sort(universal_sig_sites))

# Using escape_mat_cluster (matrix: rows = companies, cols = sites)
############## CLUSTER 13 ########################


# 1. Cluster-level: 2 SD threshold
mean_escape_per_site <- colMeans(escape_mat_cluster_13, na.rm = TRUE)
mean_escape_per_site <- mean_escape_per_site[!is.nan(mean_escape_per_site)]

overall_mean <- mean(mean_escape_per_site)
site_sd <- sd(mean_escape_per_site)
cluster_threshold <- overall_mean + 2 * site_sd

# LIST 1: Cluster-only (2SD threshold)
cluster_only_sites <- names(mean_escape_per_site)[mean_escape_per_site >= cluster_threshold]

# 2. Per-company stats (0.5 SD threshold)
per_company_stats <- lapply(1:nrow(escape_mat_cluster_13), function(i) {
  company_row <- escape_mat_cluster_13[i, , drop = FALSE]
  company_mean <- mean(company_row, na.rm = TRUE)
  company_sd <- sd(company_row, na.rm = TRUE)
  list(mean = company_mean, sd = company_sd)
})

# LIST 2: Universal sites (2SD cluster + 0.5SD EVERY company)
universal_sig_sites <- cluster_only_sites
for(i in 1:nrow(escape_mat_cluster_13)) {
  company_row <- as.numeric(escape_mat_cluster_13[i, cluster_only_sites])
  company_threshold <- per_company_stats[[i]]$mean + 0.5 * per_company_stats[[i]]$sd
  
  universal_sig_sites <- universal_sig_sites[company_row >= company_threshold]
}

# Summary
cat("=== CLUSTER 13 EPITOPE SITES ===\n")
cat("Cluster threshold (2 SD):", round(cluster_threshold, 4), "\n")
cat("\nLIST 1 - Cluster-only (2SD):", length(cluster_only_sites), "sites\n")
print(sort(cluster_only_sites))

cat("\nLIST 2 - Universal (2SD cluster + 0.5SD all companies):", 
    length(universal_sig_sites), "sites\n")
print(sort(universal_sig_sites))

# Using escape_mat_cluster (matrix: rows = companies, cols = sites)
############## CLUSTER 14 ########################


# 1. Cluster-level: 2 SD threshold
mean_escape_per_site <- colMeans(escape_mat_cluster_14, na.rm = TRUE)
mean_escape_per_site <- mean_escape_per_site[!is.nan(mean_escape_per_site)]

overall_mean <- mean(mean_escape_per_site)
site_sd <- sd(mean_escape_per_site)
cluster_threshold <- overall_mean + 2 * site_sd

# LIST 1: Cluster-only (2SD threshold)
cluster_only_sites <- names(mean_escape_per_site)[mean_escape_per_site >= cluster_threshold]

# 2. Per-company stats (0.5 SD threshold)
per_company_stats <- lapply(1:nrow(escape_mat_cluster_14), function(i) {
  company_row <- escape_mat_cluster_14[i, , drop = FALSE]
  company_mean <- mean(company_row, na.rm = TRUE)
  company_sd <- sd(company_row, na.rm = TRUE)
  list(mean = company_mean, sd = company_sd)
})

# LIST 2: Universal sites (2SD cluster + 0.5SD EVERY company)
universal_sig_sites <- cluster_only_sites
for(i in 1:nrow(escape_mat_cluster_14)) {
  company_row <- as.numeric(escape_mat_cluster_14[i, cluster_only_sites])
  company_threshold <- per_company_stats[[i]]$mean + 0.5 * per_company_stats[[i]]$sd
  
  universal_sig_sites <- universal_sig_sites[company_row >= company_threshold]
}

# Summary
cat("=== CLUSTER 14 EPITOPE SITES ===\n")
cat("Cluster threshold (2 SD):", round(cluster_threshold, 4), "\n")
cat("\nLIST 1 - Cluster-only (2SD):", length(cluster_only_sites), "sites\n")
print(sort(cluster_only_sites))

cat("\nLIST 2 - Universal (2SD cluster + 0.5SD all companies):", 
    length(universal_sig_sites), "sites\n")
print(sort(universal_sig_sites))


# Z-normalize ALL columns (sites across all antibodies) for analysis of escape scores
# in lineages
combined_num_normalized <- as.data.frame(scale(combined_num))

# Check it worked (means ~0, SDs ~1 per column)
colMeans(combined_num_normalized, na.rm=TRUE)  # Should be ~0
apply(combined_num_normalized, 2, sd, na.rm=TRUE)  # Should be ~1

#install.packages("clipr")
library(clipr)
column_names <- colnames(combined_num_normalized) 
clipr::write_clip(column_names)







