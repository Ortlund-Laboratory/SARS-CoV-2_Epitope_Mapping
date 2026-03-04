library(ggplot2)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(tidyr)
library(seqinr)
require(ggseqlogo)
library(tidyverse)
library(grid)
library(readr)
library(caret)

############################
# File I/O
############################

path <- "/home/adkeith@Eu.Emory.Edu/DMS_Workflow/nucleocapsid/Covid_Update_Correlogram/consistent_heatmaps"

files <- list.files(
  path = path,
  pattern = "_escape_fractions_raw\\.csv$",
  full.names = TRUE
)

names_vec <- sub("_escape_fractions_raw\\.csv$", "", basename(files))

dat_list <- lapply(files, function(f) {
  read.csv(f, stringsAsFactors = FALSE)
})
names(dat_list) <- names_vec

for (i in seq_along(dat_list)) {
  assign(names_vec[i], dat_list[[i]])
}

reference <- read.csv("reference.csv", header = TRUE)

############################
# Combine data
############################

df_list <- mget(names_vec)
names(df_list) <- names_vec

for (nm in names(df_list)) {
  names(df_list[[nm]])[names(df_list[[nm]]) == "Escape_Score"] <- nm
}

reference$row_id <- 1:nrow(reference)
combined <- reference

for (d in df_list) {
  combined <- merge(
    combined,
    d,
    by   = "Mutation",
    all.x = TRUE,
    all.y = FALSE,
    sort = FALSE
  )
}

combined <- combined[order(combined$row_id), ]
combined$row_id <- NULL
rownames(combined) <- NULL

combined_num <- combined[, sapply(combined, is.numeric)]

process <- preProcess(as.data.frame(combined_num), method = c("range"))
combined_num_normalized <- predict(process, as.data.frame(combined_num))

############################
# Mutation parsing helpers
############################

fun_first <- function(x) substring(x, 1, 1)
fun_last  <- function(x) substr(x, nchar(x), nchar(x))
fun_site  <- function(x) substr(x, 2, nchar(x) - 1)

tmp <- data.frame(
  wildtype = sapply(reference$Mutation, fun_first),
  site     = as.numeric(sapply(reference$Mutation, fun_site)),
  mutation = sapply(reference$Mutation, fun_last)
)
colnames(tmp) <- c("wildtype", "site", "mutation")
frames <- distinct(tmp)
frames$site <- as.integer(frames$site)

# amino acid order
polar    <- c("H", "C", "S", "T", "N", "Q")
nonpolar <- c("G", "A", "V", "L", "I", "M", "P")
aromatic <- c("F", "Y", "W")
positive <- c("K", "R")
negative <- c("D", "E")
aa_order <- c(negative, positive, polar, nonpolar, aromatic)

############################
# Heatmaps
############################

escape_trimmed <- data.frame(
  site     = tmp$site,
  mutation = tmp$mutation
)

cols <- colnames(combined_num_normalized)

out_dir <- "escape_heatmaps"
dir.create(out_dir, showWarnings = FALSE)

make_heatmaps_for_col <- function(col_name) {
  escape_trimmed$Escape <- combined_num_normalized[[col_name]]
  antibody <- col_name
  
  ## --- Residues 1–209 ---
  start <- 1
  end   <- 209
  
  mut_range <- escape_trimmed[escape_trimmed$site >= start &
                                escape_trimmed$site <= end,
                              c("site", "mutation", "Escape")]
  mut_range$site <- as.factor(mut_range$site)
  
  frames_range <- subset(frames, site >= start & site <= end)
  frames_range$site <- as.factor(frames_range$site)
  
  g1 <- ggplot(mut_range, aes(site, mutation, size = Escape * Escape * Escape)) +
    geom_tile(color = "white",
              fill  = "#FCF0F0",
              lwd   = 0.1,
              linetype = 1,
              alpha = 1) +
    ylim(rev(aa_order)) +
    geom_tile(data = frames_range,
              size   = 0,
              height = 1,
              fill   = "white",
              colour = "white") +
    geom_point(aes(colour = Escape), alpha = 1) +
    scale_colour_distiller(
      palette   = "RdPu",
      direction = +1,
      na.value  = "#FCF0F0",
      limits    = c(0, max(escape_trimmed$Escape, na.rm = TRUE)),
      name      = "Escape"
    ) +
    # remove size legend entirely
    scale_size(range = c(0, 2), guide = "none") +
    geom_point(inherit.aes = FALSE,
               data = frames_range,
               aes(site, wildtype),
               shape  = 16,
               size   = 1.25,
               colour = "black") +
    xlab("") + ylab("") +
    scale_x_discrete(breaks = as.character(seq(5, 205, by = 5))) +
    theme(
      panel.border     = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # bigger legend text and key
      legend.title = element_text(size = 10, face = "bold"),
      legend.text  = element_text(size = 9),
      legend.key.height = unit(0.35, "inch"),
      # place legend on right, centered vertically
      legend.position = c(1.02, 0.5),
      legend.justification = c("left", "center")
    )
  
  ggsave(
    filename = file.path(out_dir, paste0(antibody, "_EscapeFraction_heatmap01.png")),
    plot     = g1,
    width    = 16, height = 3, dpi = 300, units = "in", device = "png"
  )
  
  ## --- Residues 210–420 ---
  start <- 210
  end   <- 420
  
  mut_range <- escape_trimmed[escape_trimmed$site >= start &
                                escape_trimmed$site <= end,
                              c("site", "mutation", "Escape")]
  mut_range$site <- as.factor(mut_range$site)
  
  frames_range <- subset(frames, site >= start & site <= end)
  frames_range$site <- as.factor(frames_range$site)
  
  g2 <- ggplot(mut_range, aes(site, mutation, size = Escape * Escape * Escape)) +
    geom_tile(color = "white",
              fill  = "#FCF0F0",
              lwd   = 0.1,
              linetype = 1,
              alpha = 1) +
    ylim(rev(aa_order)) +
    geom_tile(data = frames_range,
              size   = 0,
              height = 1,
              fill   = "white",
              colour = "white") +
    geom_point(aes(colour = Escape), alpha = 1) +
    scale_colour_distiller(
      palette   = "RdPu",
      direction = +1,
      na.value  = "#FCF0F0",
      limits    = c(0, max(escape_trimmed$Escape, na.rm = TRUE)),
      name      = "Escape"
    ) +
    scale_size(range = c(0, 2), guide = "none") +
    geom_point(inherit.aes = FALSE,
               data = frames_range,
               aes(site, wildtype),
               shape  = 16,
               size   = 1.25,
               colour = "black") +
    xlab("") + ylab("") +
    scale_x_discrete(breaks = as.character(seq(210, 415, by = 5))) +
    theme(
      panel.border     = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text  = element_text(size = 9),
      legend.key.height = unit(0.35, "inch"),
      legend.position = c(1.02, 0.5),
      legend.justification = c("left", "center")
    )
  
  ggsave(
    filename = file.path(out_dir, paste0(antibody, "_EscapeFraction_heatmap02.png")),
    plot     = g2,
    width    = 16, height = 3, dpi = 300, units = "in", device = "png"
  )
}

for (col in cols) {
  make_heatmaps_for_col(col)
}

############################
# Per‑site escapes (no fixed blocks)
############################

# DIAGNOSTIC - run this and tell me the output
print("=== BLOCK COUNT ===")
print(nrow(avg_mat))
print("=== FIRST 10 SITES ===")
print(head(avg_mat$site_num, 10))
print("=== SITES 208-212 ===")
print(avg_mat$site_num[207:212])  # should show 208,209,210,211,212
print("=== 1-209 COUNT ===")
print(sum(avg_mat$site_num >= 1 & avg_mat$site_num <= 209, na.rm = TRUE))
print("=== HOW MANY AT EXACTLY 209? ===")
print(sum(avg_mat$site_num == 209, na.rm = TRUE))


# site for each row in combined_num_normalized
site_of_row <- as.numeric(sapply(reference$Mutation, fun_site))

# sanity: this should be same length as number of rows
stopifnot(length(site_of_row) == nrow(combined_num_normalized))

# per‑site means for each antibody column
avg_mat <- combined_num_normalized %>%
  mutate(site = site_of_row) %>%
  group_by(site) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  ungroup()

# make sure site is integer and sorted
avg_mat$site <- as.integer(avg_mat$site)
avg_mat <- avg_mat[order(avg_mat$site), ]

cols <- setdiff(colnames(avg_mat), "site")

out_dir <- "escape_barplots"
dir.create(out_dir, showWarnings = FALSE)

for (col in cols) {
  # thresholds based on per‑site means
  m_   <- mean(avg_mat[[col]], na.rm = TRUE)
  sd_  <- sd(avg_mat[[col]], na.rm = TRUE)
  cut_1 <- m_ + 1 * sd_
  cut_2 <- m_ + 2 * sd_
  
  # second criterion: count original values per site above mean+1SD
  col_idx <- which(colnames(combined_num_normalized) == col)
  
  orig_threshold <- mean(combined_num_normalized[, col_idx], na.rm = TRUE) +
    sd(combined_num_normalized[, col_idx], na.rm = TRUE)
  
  # logical: TRUE if value for that row >= orig_threshold
  above_vec <- combined_num_normalized[, col_idx] >= orig_threshold
  
  # count per site (names will be the site value)
  count_above_orig <- tapply(
    above_vec,
    site_of_row,
    function(x) sum(x, na.rm = TRUE)
  )
  
  # align counts with avg_mat$site
  count_aligned <- count_above_orig[as.character(avg_mat$site)]
  
  # per‑site escape & level
  avg_mat$escape <- ifelse(avg_mat[[col]] < 0, 0, avg_mat[[col]])
  avg_mat$level  <- ifelse(avg_mat[[col]] < cut_1, "grey",
                           ifelse(avg_mat[[col]] < cut_2, "red_1", "red_2"))
  avg_mat$level  <- as.factor(avg_mat$level)
  
  # red stars: at least 8 mutations in that site above orig_threshold AND mean > cut_1
  avg_mat$red_star <- count_aligned >= 8 & avg_mat[[col]] > cut_1 & !is.na(avg_mat[[col]])
  
  # export top escapes and star sites
  escapees_C02 <- subset(avg_mat, avg_mat[[col]] > cut_2, select = c("site", col))
  top_escapes <- na.omit(escapees_C02)
  write.csv(top_escapes,
            file = file.path(out_dir, paste0(col, "_top_escapes.csv")),
            row.names = FALSE)
  
  red_star_sites <- avg_mat[avg_mat$red_star, c("site", col)]
  red_star_sites$mean_value       <- avg_mat[[col]][avg_mat$red_star]
  red_star_sites$count_above_orig <- count_aligned[avg_mat$red_star]
  write.csv(red_star_sites,
            file = file.path(out_dir, paste0(col, "_red_star_sites.csv")),
            row.names = FALSE)
  
  #### Residues 1–209
  # because of xlim need to do 1:210
  average_escape1 <- subset(avg_mat, site >= 1 & site <= 209)
  
  p1 <- ggplot(average_escape1, aes(site, escape, fill = level)) +
    geom_bar(stat = "identity") +
    geom_segment(aes(x = 1, xend = 209, y = cut_1, yend = cut_1),
                 linetype = 3, size = 0.15) +
    geom_segment(aes(x = 1, xend = 209, y = cut_2, yend = cut_2),
                 linetype = 3, size = 0.15) +
    geom_point(data = average_escape1[average_escape1$red_star, ],
               aes(x = site, y = escape + 0.2 * max(escape, na.rm = TRUE)),
               shape = 8, size = 1, color = "red", stroke = 0.8) +
    scale_fill_manual(values = c("grey"  = "#d9d9d9",
                                 "red_1" = "#dd3497",
                                 "red_2" = "#4a2267")) +
    xlim(1, 210) +
    theme_void() +
    theme(legend.position = "none")
  
  ggsave(filename = file.path(out_dir, paste0(col, "_average_escape_01.png")),
         plot = p1, width = 12, height = 0.4, dpi = 300,
         units = "in", device = "png")
  
  #### Residues 210–420
  # because of xlim extend to 209:420
  average_escape2 <- subset(avg_mat, site >= 210 & site <= 420)
  
  p2 <- ggplot(average_escape2, aes(site, escape, fill = level)) +
    geom_bar(stat = "identity") +
    geom_segment(aes(x = 210, xend = 420, y = cut_1, yend = cut_1),
                 linetype = 3, size = 0.15) +
    geom_segment(aes(x = 210, xend = 420, y = cut_2, yend = cut_2),
                 linetype = 3, size = 0.15) +
    geom_point(data = average_escape2[average_escape2$red_star, ],
               aes(x = site, y = escape + 0.2 * max(escape, na.rm = TRUE)),
               shape = 8, size = 1, color = "red", stroke = 0.8) +
    scale_fill_manual(values = c("grey"  = "#d9d9d9",
                                 "red_1" = "#dd3497",
                                 "red_2" = "#4a2267")) +
    xlim(209, 420) +
    theme_void() +
    theme(legend.position = "none")
  
  ggsave(filename = file.path(out_dir, paste0(col, "_average_escape_02.png")),
         plot = p2, width = 12, height = 0.4, dpi = 300,
         units = "in", device = "png")
  
  # clean temporary columns before next antibody
  avg_mat$escape   <- NULL
  avg_mat$level    <- NULL
  avg_mat$red_star <- NULL
}
