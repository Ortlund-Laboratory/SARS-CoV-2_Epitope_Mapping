library(ggplot2)
#library(hrbrthemes)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(seqinr)
require(ggseqlogo)
library(tidyverse)
library(grid)
library(readr)

################################################################################
###Rename 
################################################################################
antibody <- "Ab339"
ab_name <- antibody

escape_counts <- read_tsv("/home/adkeith@Eu.Emory.Edu/DMS_Workflow/nucleocapsid/test_p23068_analysis_Ab339_CorDx/results/escape_scores/esc_variant_counts.txt")
input_counts <- read_tsv("/home/adkeith@Eu.Emory.Edu/DMS_Workflow/nucleocapsid/test_p23068_analysis_Ab339_CorDx/results/escape_scores/ref_variant_counts.txt")

escape_counts <- escape_counts[,2:3]
input_counts <- input_counts[,2:3]

################################################################################
### Histogram for read counts
################################################################################
ggsave(filename = paste(antibody, "_Number_of_Reads_escaped.png", sep=""), 
       ggplot(escape_counts, aes(x=count), title = "Mutations reads")+
         geom_histogram(color = "black", fill = "grey", ) +
         xlim(0,quantile(escape_counts$count, probs = .95, na.rm = T)) +
         ylim(0,5000)+
         # geom_vline(aes(xintercept=mu_bc_depth),
         #            linetype="dashed",
         #            color = "red") + 
         #geom_text(aes(label = mu), x = 3, y = 3, vjust = "inward", hjust = "inward") +
         ggtitle("Escape Read Counts") +
         xlab("# Reads") + 
         ylab("Count") +
         theme_bw(base_size = 10),
       width = 3, height = 2, dpi = 300, units = "in", device='png')

ggsave(filename = paste(antibody, "_Number_of_Reads_Myc.png", sep=""), 
       ggplot(input_counts, aes(x=count), title = "Mutations reads")+
         geom_histogram(color = "black", fill = "grey") +
         xlim(0,quantile(input_counts$count, probs = .99, na.rm = T)) +
         ylim(0,30000)+
         # geom_vline(aes(xintercept=mu_bc_depth),
         #            linetype="dashed",
         #            color = "red") + 
         #geom_text(aes(label = mu), x = 3, y = 3, vjust = "inward", hjust = "inward") +
         ggtitle("Reference Read Counts") +
         xlab("# Reads") + 
         ylab("Count") +
         theme_bw(base_size = 10),
       width = 3, height = 2, dpi = 300, units = "in", device='png')

#####################################################################
#Combine duplicated mutations (same mutation with multiple barcodes)
#####################################################################
escape_combined <- escape_counts %>% 
  group_by(mutation) %>% 
  summarise_all(funs(sum))

reference_combined <- input_counts %>% 
  group_by(mutation) %>% 
  summarise_all(funs(sum))

colnames(reference_combined) <- c("name", "depth")
colnames(escape_combined) <- c("name", "depth")

####################################################################
##Renaming for escape
####################################################################

fun_first <- function(x) {
  substring(x, 1,1)
}

fun_last <- function(x){
  substr(x, nchar(x), nchar(x))
}

fun_site <- function(x){
  substr(x, 2, nchar(x)-1)
}

wildtype <- data.frame(sapply(escape_combined[1], fun_first))
#colnames(wildtype) <- c("wildtype")
escape_combined$wildtype <- wildtype[[1]]

mutation <- data.frame(sapply(escape_combined[1], fun_last))
colnames(mutation) <- c("mutation")
escape_combined$mutation <- as.factor(mutation[[1]])

site <- data.frame(sapply(escape_combined[1], fun_site))
colnames(site) <- c("site")
escape_combined$site <- as.numeric(site[[1]])
#Reorder the columns (name, wild type, site, mutation, depth)
escape_combined <- escape_combined[,c(1,3,5,4,2)]

###Order reference_combined by site and mutation
escape_combined <- escape_combined[
  with(escape_combined, order(site, mutation)),
]



####################################################################
##Renaming for reference
####################################################################
wildtype <- data.frame(sapply(reference_combined[1], fun_first))
#colnames(wildtype) <- c("wildtype")
reference_combined$wildtype <- wildtype[[1]]

mutation <- data.frame(sapply(reference_combined[1], fun_last))
colnames(mutation) <- c("mutation")
reference_combined$mutation <- as.factor(mutation[[1]])

site <- data.frame(sapply(reference_combined[1], fun_site))
colnames(site) <- c("site")
reference_combined$site <- as.numeric(site[[1]])
#Reorder the columns (name, wild type, site, mutation, depth)
reference_combined <- reference_combined[,c(1,3,5,4,2)]

###Order reference_combined by site and mutation
reference_combined <- reference_combined[
  with(reference_combined, order(site, mutation)),
]

#Clean-up
remove(mutation)
remove(site)
remove(wildtype)




###Remove stop codons
escape_trimmed <- escape_combined[which(escape_combined$mutation != "*"),]
reference_trimmed <- reference_combined[which(reference_combined$mutation != "*"),]


###############################################################################
###Add missing residues
###############################################################################
Wuhan <- read.fasta("N_Wuhan.fasta")

###First reference data
### If there are residues without any data, these have to be added manually
#1. Define the range of amino acids present
seq_range <- min(reference_trimmed$site):max(reference_trimmed$site)
missing_aa <- seq_range[!seq_range %in% unique(reference_trimmed$site)]
#2. Add the missing residue(s)
missing_data <- data.frame(site = missing_aa, 
                           name = paste(toupper(Wuhan[["N_Wuhan"]][missing_aa]), missing_aa, toupper(Wuhan[["N_Wuhan"]][missing_aa])),
                           mutation = toupper(Wuhan[["N_Wuhan"]][missing_aa]), 
                           wildtype = toupper(Wuhan[["N_Wuhan"]][missing_aa]))
complete_reference <- rbind.fill(reference_trimmed, missing_data)

#Expand the dataset to include NA values for synonymous amino acids
all <- complete_reference %>% expand(site, mutation)
# join with all, n will be NA for obs. in all that are not present in v
reference_trimmed = complete_reference %>% group_by_at(vars(wildtype, site, mutation)) %>% 
  right_join(all)


###Next escape data
### If there are residues without any data, these have to be added manually
#1. Define the range of amino acids present
seq_range <- min(escape_trimmed$site):max(escape_trimmed$site)
missing_aa <- seq_range[!seq_range %in% unique(escape_trimmed$site)]
#2. Add the missing residue(s)
missing_data <- data.frame(site = missing_aa, 
                           name = paste(toupper(Wuhan[["N_Wuhan"]][missing_aa]), missing_aa, toupper(Wuhan[["N_Wuhan"]][missing_aa])),
                           mutation = toupper(Wuhan[["N_Wuhan"]][missing_aa]), 
                           wildtype = toupper(Wuhan[["N_Wuhan"]][missing_aa]))
complete_escape <- rbind.fill(escape_trimmed, missing_data)

#Expand the dataset to include NA values for synonymous amino acids
all <- complete_escape %>% expand(site, mutation)
# join with all, n will be NA for obs. in all that are not present in v
escape_trimmed = complete_escape %>% group_by_at(vars(wildtype, site, mutation)) %>% 
  right_join(all)



###Order reference_trimmed and escape_trimmed by site and mutation
escape_trimmed <- escape_trimmed[
  with(escape_trimmed, order(site, mutation)),
]
reference_trimmed <- reference_trimmed[
  with(reference_trimmed, order(site, mutation)),
]


###################################################################
### Calculate abundance: n/N
###################################################################
escape_trimmed$abundance <- escape_trimmed$depth/sum(escape_trimmed$depth, na.rm=TRUE)
reference_trimmed$abundance <- reference_trimmed$depth/sum(reference_trimmed$depth, na.rm=TRUE)

###Then adjust the lowest reference abundance values
#Use 95 percentile
#This helps remove exaggerated escape fractions caused by dividing by a very small number
cutoff  <- quantile(reference_trimmed$abundance, probs = .05, na.rm = T)
reference_trimmed$abundance <- ifelse(reference_trimmed$abundance<cutoff, cutoff, reference_trimmed$abundance)

##Fishers exact test
esc_fisher <- escape_trimmed[,1:5]
esc_fisher$ref_depth <- reference_trimmed$depth
esc_sum <- sum(na.omit(esc_fisher$depth))
ref_sum <- sum(na.omit(esc_fisher$ref_depth))

p_values <- data.frame(matrix(ncol = 1, nrow = nrow(esc_fisher)))
colnames(p_values) <- "p"
#Calculate p values for each comparison (using Fisher's exact test)
for (i in 1:nrow(esc_fisher)){
  if (!anyNA(esc_fisher[i,])){
    f_ <- matrix(unlist(c(esc_fisher[i,5], 
                          esc_fisher[i,6],
                          esc_sum-esc_fisher[i,5],
                          ref_sum-esc_fisher[i,6])),2)
    p_ <- fisher.test(f_)
    p_values$p[i] <-p_$p.value
  }
}

p_values$p_adj <- p.adjust(p_values$p, "bonferroni")
p_values$p_adj <- ifelse(p_values$p_adj == 0, 1e-300, p_values$p_adj)

###################################################################
### Calculate escape fractions
### And add adjusted p values
###################################################################
escape_trimmed$Escape <- (escape_trimmed$abundance/reference_trimmed$abundance)
escape_trimmed$p_adj <- p_values$p_adj
max(escape_trimmed$Escape, na.rm=T)
min(escape_trimmed$Escape, na.rm=T)
#Write out raw escape score before normalizing
escape_raw_out <- na.omit(escape_trimmed[,c(1, ncol(escape_trimmed)-2)])
colnames(escape_raw_out) <- c("Mutation", "EscapeScore")
#write.csv(escape_raw_out, file = paste(antibody,"_raw_escape.csv", sep = ""), row.names = FALSE)

escape_trimmed$Non_Normalized_Enrich <- escape_trimmed$Escape

###Normalize the data between the 99 and 1 percentiles
max(escape_trimmed$Escape, na.rm=T)
upper_limit  <- quantile(escape_trimmed$Escape, probs = .99, na.rm = T)
lower_limit <- min(escape_trimmed$Escape, na.rm=T)
escape_trimmed$Escape <- ifelse(escape_trimmed$Escape>upper_limit, upper_limit, escape_trimmed$Escape)
escape_trimmed$Escape <- (escape_trimmed$Escape-lower_limit)/(upper_limit-lower_limit)

ggsave(filename = paste(antibody,"_EscapeHistogram_preArcSine.png", sep=""), 
       ggplot(escape_trimmed, aes(x=Escape), title = "Escape Fractions")+
         #geom_histogram(aes(y=..density..), color = "black", fill = "grey", binwidth = 0.025) +
         geom_histogram(color = "black", fill = "grey") +
         xlim(0,quantile(escape_trimmed$Escape, probs = .99, na.rm = T)) +
         xlab("Escape Fraction") + 
         ylab("Count") +
         theme_bw(base_size = 10),
       width = 3, height = 2, dpi = 300, units = "in", device='png')


################################################################################
### Generate a matrix for plotting sequence logos
### This matrix is also used to calculate average escape scores
################################################################################
#remove the extra rwos that only contain the site number and nothing else 
escape_trimmed <- escape_trimmed[which(!is.na(escape_trimmed$mutation)),]

logo_matrix <- matrix(ncol = nrow(escape_trimmed)/20,nrow=20)
row.names(logo_matrix) <- escape_trimmed$mutation[1:20]
colnames(logo_matrix) <- seq(2,ncol(logo_matrix)+1, 1)


for(i in 0:ncol(logo_matrix)-1){
  logo_matrix[,i+1] <- escape_trimmed$Escape[seq(from = 20*i+1, to=20*i+20, by = 1)]
  
  ################################################################################
  ### Generate output files for per-variant and per-site escapes
  ################################################################################
  average_escape <- as.data.frame(colMeans(logo_matrix, na.rm=T))
  average_escape$site <- seq(2,419,1)
  colnames(average_escape) <- c("escape", "site")
  
  #Individual escapes:
  escape_fractions_out <- na.omit(escape_trimmed[,c(1,7,8)])
  #escape_fractions_out$name <- paste(escape_fractions_out$wildtype, escape_fractions_out$site, escape_fractions_out$mutation, sep = "")
  write.csv(escape_fractions_out, file = paste(antibody,"_escape_fractions.csv", sep = ""), row.names = FALSE)
  
  #Average escapes:
  average_escape_out <- as.data.frame(average_escape$site)
  average_escape_out$'Average Escape Score' <- average_escape$escape
  colnames(average_escape_out) <- c("Site", "Average Escape Score")
  write.csv(average_escape_out, file = paste(antibody,"_average_escape.csv", sep=""), row.names = FALSE)
  
}
