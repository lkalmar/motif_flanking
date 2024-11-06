library(stringr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggseqlogo)

#########
### Masking process (optional, fully masked data table is available to load in)
#########

data <- read.csv("motif_data.csv")
masking_table <- read.csv("masking_table.csv", row.names = 1)
data$Masked_motif <- NA
for (i in 1:length(rownames(data))){
  if (!is.na(masking_table[data$up_to_date_class[i], "masking_regex"])){
    data$Masked_motif[i] <- gsub(masking_table[data$up_to_date_class[i], "masking_regex"],
                              masking_table[data$up_to_date_class[i], "replacement_string"],
                              data$Motif[i])
  }
}
write.csv(data, "motif_data.csv", row.names = FALSE)

#########
### Manual masking of missing masked motifs
#########


#########
### Read in the finished masked motif data instead
#########

data <- read.csv("motif_data_fixed.csv")
sum(is.na(data$Masked_motif))
data$Motif_variable_content <- gsub("_", "", data$Masked_motif)

aaindex_table <- read.table("aaindex.txt", sep = "*", header = TRUE, row.names = 1)

calc_mean <- function(indexAC, seq){
  vec <- unlist(strsplit(seq, ""))
  sum <- 0
  for (i in vec){
    sum <- sum + as.numeric(aaindex_table[indexAC, i])
  }
  return(sum / length(vec))
}

#########
### Calculating mean values for different regions
#########

minlen <- 3
data <- data[,1:14]
for (aaID in rownames(aaindex_table)){
  data[, paste0(aaID,"_N11to20")] <- NA
  data[, paste0(aaID,"_N1to10")] <- NA
  data[, paste0(aaID,"_N1to5")] <- NA
  data[, paste0(aaID,"_MotifOrig")] <- NA
  data[, paste0(aaID,"_MotifMasked")] <- NA
  data[, paste0(aaID,"_C11to20")] <- NA
  data[, paste0(aaID,"_C1to10")] <- NA
  data[, paste0(aaID,"_C1to5")] <- NA
  for (i in 1:length(rownames(data))){
    # N terminal
    if (nchar(data$Flanking_N[i]) == 20){
      data[i,paste0(aaID,"_N11to20")] <- calc_mean(aaID, substr(data$Flanking_N[i], 1, nchar(data$Flanking_N[i])-10))
    }
    if (nchar(data$Flanking_N[i]) >= 10){
      start <- nchar(data$Flanking_N[i]) - 9
      if (start < 1){
        start <- 1
      }
      data[i,paste0(aaID,"_N1to10")] <- calc_mean(aaID, substr(data$Flanking_N[i], start, nchar(data$Flanking_N[i])))
    }
    if (nchar(data$Flanking_N[i]) >= 5){
      start <- nchar(data$Flanking_N[i]) - 4
      if (start < 1){
        start <- 1
      }
      data[i,paste0(aaID,"_N1to5")] <- calc_mean(aaID, substr(data$Flanking_N[i], start, nchar(data$Flanking_N[i])))
    }
    # Motif
    if (nchar(data$Motif[i]) >= minlen){
      data[i,paste0(aaID,"_MotifOrig")] <- calc_mean(aaID, data$Motif[i])
    }
    if (nchar(data$Motif_variable_content[i]) >= minlen){
      data[i,paste0(aaID,"_MotifMasked")] <- calc_mean(aaID, data$Motif_variable_content[i])
    }
    #C terminal
    if (nchar(data$Flanking_C[i]) == 20){
      data[i,paste0(aaID,"_C11to20")] <- calc_mean(aaID, substr(data$Flanking_C[i], 11, nchar(data$Flanking_C[i])))
    }
    if (nchar(data$Flanking_C[i]) >= 10){
      end <- 10
      data[i,paste0(aaID,"_C1to10")] <- calc_mean(aaID, substr(data$Flanking_C[i], 1, end))
    }
    if (nchar(data$Flanking_C[i]) >= 5){
      end <- 5
      data[i,paste0(aaID,"_C1to5")] <- calc_mean(aaID, substr(data$Flanking_C[i], 1, end))
    }
  }
}



write.csv(data, "motif_data.csv", row.names = FALSE)
data <- read.csv("motif_data.csv")

#########
### Get disordered regions from MobiDB (minimum length 10aa)
#########

extract_segments <- function(mask, sequence, min_length = 10) {
  mask_vec <- unlist(strsplit(mask, ""))
  sequence_vec <- unlist(strsplit(sequence, ""))
  segments <- list()
  rle_mask <- rle(mask_vec)
  start_positions <- cumsum(c(1, head(rle_mask$lengths, -1)))
  ones_positions <- start_positions[rle_mask$values == "1"]
  ones_lengths <- rle_mask$lengths[rle_mask$values == "1"]
  for (i in seq_along(ones_positions)) {
    if (ones_lengths[i] >= min_length) {
      start <- ones_positions[i]
      end <- start + ones_lengths[i] - 1
      segments[[length(segments) + 1]] <- paste(sequence_vec[start:end], collapse = "")
    }
  }
  return(segments)
}

mobiDBdata <- readLines("mobidb_search_2024-10-26T15-52-49.fasta")
disordered_regions <- c()
pass <- 0
for (line in mobiDBdata){
  if (pass == 1){
    seq <- line
  }
  if (pass == 2){
    disordered_regions <- c(disordered_regions,
                            unlist(extract_segments(line,seq)))
  }
  pass <- 0
  if (grepl("\\|sequence", line)){
    pass <- 1
  }
  if (grepl("\\|curated-disorder-merge", line)){
    pass <- 2
  }
}

#########
### Randomly select 1000 regions, extract random 10-mers, calculate aaindex means
#########

set.seed(123)
selected_fragments <- sample(disordered_regions, 1000, replace = FALSE)
extracted_segments <- sapply(selected_fragments, function(fragment) {
  if (nchar(fragment) >= 10) {
    start_pos <- sample(1:(nchar(fragment) - 9), 1)
    substring(fragment, start_pos, start_pos + 9)
  } else {
    NA
  }
})
extracted_segments <- extracted_segments[!is.na(extracted_segments)]

disordered_df <- data.frame(Sequence = extracted_segments,
                            Source = "MobiDB")
for (aaID in rownames(aaindex_table)){
  disordered_df[, paste0(aaID,"_MobiDB")] <- NA
  for (i in 1:length(rownames(disordered_df))){
    disordered_df[i, paste0(aaID,"_MobiDB")] <- calc_mean(aaID, disordered_df$Sequence[i])
  }
}

#########
### Creating long format tables for visualisation (multiple options below to create visdf)
#########

### Full data set
visdf <- melt(data[,15:78])
visdf <- rbind(visdf,
               melt(disordered_df[,3:10]))
visdf$AAindex <- str_split_i(visdf$variable, "_", 1)
visdf$Region <- str_split_i(visdf$variable, "_", 2)
visdf <- visdf[!is.na(visdf$value),]

### Filtered data set to have one instance / class
data_nr <- data[0,]
set.seed(123)
for (class in unique(data$up_to_date_class)){
  tmpdf <- subset(data, up_to_date_class == class)
  data_nr <- rbind(data_nr,
                   tmpdf[sample(length(rownames(tmpdf)),1),])
}
visdf <- melt(data_nr[,15:78])
visdf <- rbind(visdf,
               melt(disordered_df[,3:10]))
visdf$AAindex <- str_split_i(visdf$variable, "_", 1)
visdf$Region <- str_split_i(visdf$variable, "_", 2)
visdf <- visdf[!is.na(visdf$value),]

#########
### Variables to define labels and colours
#########

my_comparisons <- list()
my_comparisons[[1]] <- list(c("MobiDB", "N11to20"),
                       c("MobiDB", "N1to10"),
                       c("MobiDB", "N1to5"),
                       c("MobiDB", "C1to5"),
                       c("MobiDB", "C1to10"),
                       c("MobiDB", "C11to20"))
my_comparisons[[2]] <- list(c("MotifMasked", "C1to5"),
                        c("MotifMasked", "N1to5"),
                        c("MotifMasked", "C1to10"),
                        c("MotifMasked", "N1to10"),
                        c("MotifMasked", "C11to20"),
                        c("MotifMasked", "N11to20"))

my_limits <- c("MobiDB",
               "N11to20",
               "N1to10",
               "N1to5",
               "MotifOrig",
               "MotifMasked",
               "C1to5",
               "C1to10",
               "C11to20")
my_labels <- c("MobiDB" = "MobiBD",
               "N11to20" = "N-term\n11-20",
               "N1to10" = "N-term\n1-10",
               "N1to5" = "N-term\n1-5",
               "MotifOrig" = "Original\nmotif",
               "MotifMasked" = "Masked\nmotif",
               "C1to5" = "C-term\n1-5",
               "C1to10" = "C-term\n1-10",
               "C11to20" = "C-term\n11-20")
my_cols <- c("MobiDB" = "seagreen",
             "N11to20" = "salmon1",
             "N1to10" = "salmon2",
             "N1to5" = "salmon3",
             "MotifOrig" = "khaki3",
             "MotifMasked" = "khaki3",
             "C1to5" = "skyblue3",
             "C1to10" = "skyblue2",
             "C11to20" = "skyblue1")

#########
### For loop to create and save violoin plots for general signature detection
#########

for (akt_index in rownames(aaindex_table)){
  for (i in 1:2){
    p <- ggplot(subset(visdf, AAindex == akt_index), aes(Region, value, fill = Region)) +
      geom_violin(adjust = 1.2, draw_quantiles = c(0.25,0.5,0.75)) +
      scale_x_discrete(limits = my_limits,
                       labels = my_labels) +
      theme_bw() +
      theme(axis.text.y = element_text(size = 11),
            axis.text.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold"),
            axis.title.x = element_blank(),
            title = element_text(size = 16, face = "bold"),
            legend.position = "none") +
      ylab("Mean of simplified AA index values") +
      ggtitle(aaindex_table[akt_index,"Description"]) + 
      scale_fill_manual(values = my_cols) +
      stat_compare_means(comparisons = my_comparisons[[i]])
    ggsave(paste0("plots/Fig3_",akt_index,"_comp_",i,".jpg"),p, width =7, height = 7)
  }
}

#########
### Intra-class comparisons and visualisations
#########
hist(unname(table(data$up_to_date_class)), breaks = 50)

#########
### Filter classes with at least 10 instances
#########
abundant_classes <- names(table(data$up_to_date_class)[table(data$up_to_date_class) >= 10])
data_abundant <- subset(data, up_to_date_class %in% abundant_classes)
visdf <- melt(data_abundant[,c(11,15:78)])
colnames(visdf)[1] <- "Source"
visdf <- rbind(visdf,
               melt(disordered_df[,2:10]))
visdf$AAindex <- str_split_i(visdf$variable, "_", 1)
visdf$Region <- str_split_i(visdf$variable, "_", 2)
visdf <- visdf[!is.na(visdf$value),]

#########
### For loop to create visualisations for class-level signature detection
#########
for (class in abundant_classes){
  fig_list <- list()
  fig_num <- 1
  akt_comparisons <- my_comparisons[[1]]
  if (grepl("LIG_PDZ", class) || grepl("LIG_WRPW", class)){
    akt_comparisons <- my_comparisons[[1]][1:3]
  }
  for (akt_index in rownames(aaindex_table)){
    p <- ggplot(subset(visdf, (AAindex == akt_index & (Source == class | Source == "MobiDB"))), aes(Region, value, fill = Region)) +
      geom_violin(adjust = 1.2, draw_quantiles = c(0.25,0.5,0.75)) +
      geom_jitter(data = subset(visdf, (AAindex == akt_index & Source == class)),alpha = 0.3, colour = "grey20", height = 0, width = 0.2) +
      scale_x_discrete(limits = my_limits,
                       labels = my_labels,
                       drop = FALSE) +
      theme_bw() +
      theme(axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 10, face = "bold"),
            axis.title.y = element_text(size = 10, face = "bold"),
            axis.title.x = element_blank(),
            title = element_text(size = 12, face = "bold"),
            legend.position = "none") +
      ylab("Mean of simplified AA index values") +
      ggtitle(aaindex_table[akt_index,"Description"]) + 
      scale_fill_manual(values = my_cols) +
      stat_compare_means(comparisons = akt_comparisons)
    fig_list[[fig_num]] <- p
    fig_num <- fig_num+1
  }
  final_plot <- ggarrange(plotlist =  fig_list, ncol = 2, nrow = 4)
  final_plot_ann <- annotate_figure(final_plot,
                                    top = text_grob(class, face = "bold", size = 18))
  ggsave(paste0("plots/class_specific/Fig_",class,"_comp_1.jpg"),final_plot_ann, width = 12, height = 15, units = "in")
}

#########
### Sequence logo visualisations for abundant classes
#########


for (class in abundant_classes){
  nterm_only <- 0
  if (grepl("LIG_PDZ", class) || grepl("LIG_WRPW", class)){
    nterm_only <- 1
  }
  # n-term
  seqvec <- data[data$up_to_date_class == class, "Flanking_N"]
  for (i in 1:length(seqvec)){
    if (nchar(seqvec[i]) < 10){
      seqvec[i] <- paste0(strrep(" ", 10 - nchar(seqvec[i])), seqvec[i])
    }
    seqvec[i] <- str_sub(seqvec[i], nchar(seqvec[i])-9, nchar(seqvec[i]))
  }
  p1 <- ggplot() + geom_logo( seqvec, method = "prob" ) + theme_logo() + theme(legend.position = "none") + ggtitle("N-terminal flanking")
  #cterm
  p2 <- ggplot() + theme_logo()
  if (nterm_only == 0){
    seqvec <- data[data$up_to_date_class == class, "Flanking_C"]
    for (i in 1:length(seqvec)){
      if (nchar(seqvec[i]) < 10){
        seqvec[i] <- paste0(seqvec[i], strrep(" ", 10 - nchar(seqvec[i])))
      }
      seqvec[i] <- str_sub(seqvec[i], 1, 10)
    }
    p2 <- ggplot() + geom_logo( seqvec, method = "prob" ) + theme_logo() + theme(legend.position = "none") + ggtitle("C-terminal flanking")
  }
  comp_plot <- ggarrange(p1, p2, ncol = 2, nrow = 1)
  final_plot <- annotate_figure(comp_plot,
                                top = text_grob(class, face = "bold", size = 18))
  ggsave(paste0("plots/seq_logos/Fig_",class,"_logo.jpg"),final_plot, width = 6, height = 3, units = "in")
}

#########
### Amino acid composition difference calculation and visualisation
#########

volume_order <- c("G","A","S","P","D","C","N","T","E","V","Q","H","M","L","I","K","R","F","Y","W")
volume_cat <- unlist(as.vector(aaindex_table["GRAR740103",2:21])[volume_order])
aa_comp <- function(seqvec, ordering){
  merged_seq <- paste0(seqvec, collapse = "")
  aacnt_vec <- unlist(strsplit(merged_seq, ""))
  abund_vec <- as.vector(unname(table(aacnt_vec)[ordering]))
  return(abund_vec/sum(abund_vec))
}

aa_comp_df <- data.frame("AminoAcid" = volume_order,
                         "Flanking_N" = aa_comp(data_nr$Flanking_N, volume_order),
                         "Flanking_C" = aa_comp(data_nr$Flanking_C, volume_order),
                         "Motif" = aa_comp(data_nr$Motif, volume_order),
                         "Motif_variable_content" = aa_comp(data_nr$Motif_variable_content, volume_order),
                         "MobiDB" = aa_comp(disordered_df$Sequence, volume_order))

for (region in colnames(aa_comp_df)[2:5]){
  diffcol <- paste0(region,"_vs_MobiDB")
  aa_comp_df[,diffcol] <- log2(aa_comp_df[,region]) - log2(aa_comp_df$MobiDB)
}

#########
### Modify aesthetics, xlab and output file name for different graphs
#########
p <- ggplot(aa_comp_df, aes(Motif_vs_MobiDB, AminoAcid, fill = as.factor(unname(volume_cat)))) +
  geom_col(colour = "black") +
  scale_y_discrete(limits = volume_order) +
  scale_fill_manual("Simplified volume category", values = c("gray80", "gray70", "gray60", "gray50", "gray40")) +
  theme_bw() +
  xlab("Difference in relative composition\nMotif vs. MobiDB (log2FC)") +
  ylab("Amino Acid (in the order of side-chain volume)") +
  theme(legend.position = "top",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_legend(nrow = 1, title.position = "top", title.hjust = 0.5))

ggsave("plots/aacont_motif_vs_mobidb.jpg", p, height = 7, width = 5)  
