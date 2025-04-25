#################################################################
###################CpG O/E among all sequences###################
#################################################################
library(seqinr)
library(ggplot2)
setwd('/Users/flickazhang/Desktop/Biohackthon_Project')

#read all sequence
for(num in 1:16){
  # Read each FASTA file
  seq_data <- read.fasta(paste0(num, "_1_sequence", ".fasta"), 
                         seqtype = "DNA", 
                         as.string = TRUE)
  seq_data <- read.fasta(paste0(num, "_2_sequence", ".fasta"), 
                         seqtype = "DNA", 
                         as.string = TRUE)
  # Create variable with dynamic name
  assign(paste0("virus_seq_20", num), seq_data)
}

#Create the Analysis Function
analyze_CpG_OE <- function(virus_seq) {
  ## Extract and prepare the sequence
  sequence <- toupper(unlist(virus_seq))
  sequence <- gsub("[^ATGC]", "", as.character(sequence))
  sequence_str <- paste(sequence, collapse = "")
  ## Find all dinucleotide repeats
  find_repeats <- function(seq_str, repeat_length) {
    motifs <- sapply(1:(nchar(seq_str) - repeat_length + 1), function(i) {
      substr(seq_str, i, i + repeat_length - 1)
    })
    repeats <- table(motifs)
    repeats[repeats > 1]  # Only keep motifs that repeat
  }
  
  di_repeats <- find_repeats(sequence_str, 2)
  di_df <- data.frame(Motif = names(di_repeats), Frequency = as.numeric(di_repeats))
  di_df <- di_df[order(-di_df$Frequency), ]
  
  ## Calculate nucleotide frequencies
  nucleotides <- table(unlist(strsplit(sequence_str, "")))
  C_freq <- nucleotides["C"] / sum(nucleotides)
  G_freq <- nucleotides["G"] / sum(nucleotides)
  
  ## Get CG frequency (handle case where CG might not exist)
  rhino_CG <- ifelse("CG" %in% di_df$Motif, 
                     di_df$Frequency[di_df$Motif == "CG"], 
                     0)
  ## Calculate O/E ratio
  OE_CpG <- (rhino_CG / sum(di_df$Frequency)) / (C_freq * G_freq)

  ## Return comprehensive results
  return(list(
    OE_ratio = OE_CpG,
    CG_count = rhino_CG,
    C_freq = C_freq,
    G_freq = G_freq,
    dinucleotide_df = di_df,
    sequence_length = length(sequence)
  ))
}

## Apply the function to all sequence
results <- list()
for(num in 1:16) {
  seq_name <- paste0("virus_seq_", num)
  results[[seq_name]] <- analyze_CpG_OE(get(seq_name))
}

#extract and compare CpG O/E
# Create summary data frame
summary_df <- data.frame(
  Sequence = names(results),
  OE_ratio = sapply(results, function(x) x$OE_ratio),
  CG_count = sapply(results, function(x) x$CG_count),
  C_freq = sapply(results, function(x) x$C_freq),
  G_freq = sapply(results, function(x) x$G_freq),
  Sequence_length = sapply(results, function(x) x$sequence_length)
)


# Sort by O/E ratio
summary_df <- summary_df[order(summary_df$OE_ratio), ]
# Vector of actual virus names in the correct order (corresponding to 1 to 16)
virus_names <- c(
  "Human rhinovirus C strain WA823M02",
  "Human rhinovirus C strain HRV-C43_p1154_sR1124_2009 polyprotein gene",
  "Rhinovirus C strain SCH-120 polyprotein gene",
  "Rhinovirus C40 isolate 470389",
  "Rhinovirus C6",
  "Human rhinovirus Ce",  # main strand
  "Human rhinovirus C strain HRV-C49_p1102_sR889_2008 polyprotein gene",
  "Rhinovirus C strain NW1",
  "Rhinovirus C isolate 1515-MY-10",
  "Rhinovirus C isolate 3430-MY-10",
  "Rhinovirus C strain 20692_95_HRV-C polyprotein gene",
  "Rhinovirus C strain 20692_39_HRV-C polyprotein gene",
  "Rhinovirus C47 isolate CA-RGDS-1001",
  "Rhinovirus C43 strain RvC43/USA/2021/WYX8S9",
  "Rhinovirus C/human/USA/MA-Broad_MGH-16204/2023 polyprotein gene",
  "Rhinovirus C/human/USA/MA-Broad_MGH-21780/2024 polyprotein gene"
)

# Sort the results list by O/E ratio
sorted_indices <- order(sapply(results, function(x) x$OE_ratio))
sorted_results <- results[sorted_indices]

# Create the summary dataframe
summary_df <- data.frame(
  Sequence = virus_names[sorted_indices],
  OE_ratio = sapply(sorted_results, function(x) x$OE_ratio),
  CG_count = sapply(sorted_results, function(x) x$CG_count),
  C_freq = sapply(sorted_results, function(x) x$C_freq),
  G_freq = sapply(sorted_results, function(x) x$G_freq),
  Sequence_length = sapply(sorted_results, function(x) x$sequence_length)
)


# Visual comparison_1
library(ggplot2)

# Optional: shorten long sequence names for better plotting (wrap at ~40 chars)
summary_df$Sequence <- stringr::str_wrap(summary_df$Sequence, width = 40)

ggplot(summary_df, aes(x = reorder(Sequence, OE_ratio), y = OE_ratio)) +
  geom_bar(stat = "identity", fill = "#66c2a5", color = "black", width = 0.7) +
  geom_text(aes(label = round(OE_ratio, 3)), 
            vjust = -0.5,            
            color = "black",         
            size = 2.5) +            
  labs(title = "CpG O/E Ratios Across Virus Sequences",
       y = "CpG O/E Ratio",
       x = "Virus Sequence") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_hline(yintercept = 0.6, linetype = "dashed", color =  "#fc8d62", size = 0.5) +
  annotate("text", x = 2, y = 0.65, label = "Strong suppression threshold", 
           color ="#fc8d62", size = 3.5) +
  ylim(0, max(summary_df$OE_ratio) * 1.1)

ggsave("CpG_OE_Ratios_Barplot.tiff", 
       width = 20, height = 6, dpi = 300, 
       bg = "white")

# Visual comparison_2
ggplot(summary_df, aes(x = reorder(Sequence, OE_ratio), y = OE_ratio)) +
  geom_bar(stat = "identity", fill = "#66c2a5", color = "black", width = 0.7) +
  # Add text labels with rounded values
  geom_text(aes(label = round(OE_ratio, 3)), 
            vjust = -0.5,            # Position above bars
            color = "black",         # Text color
            size = 2.5) +            # Text size
  labs(title = "CpG O/E Ratios Across Virus Sequences",
       y = "CpG O/E Ratio",
       x = "Virus Sequence") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # Add reference line and annotation
  geom_hline(yintercept = 0.6, linetype = "dashed", color =  "#fc8d62", size = 0.5) +
  annotate("text", x = 2, y = 0.65, label = "Strong suppression threshold", 
           color ="#fc8d62", size = 3.5) +
  # Adjust y-axis limits to accommodate text labels
  ylim(0, max(summary_df$OE_ratio) * 1.1)

#################################################################
################### CpG rate fluctuation among all sequences ###################
#################################################################
setwd('/Users/flickazhang/Desktop/Biohackthon_Project/USA_CpG_2016_2024')
fasta_files <- list.files(pattern = "_sequence\\.fasta$")
virus_sequences <- list()

for (file in fasta_files) {
  # Read the sequence
  seq_data <- read.fasta(file, seqtype = "DNA", as.string = TRUE)
  
  # Use the filename (without .fasta) as the list name
  name <- tools::file_path_sans_ext(file)
  
  # Store it in the list
  virus_sequences[[name]] <- seq_data
}

analyze_CpG_OE <- function(sequence_data) {
  ## Extract and prepare the sequence
  sequence <- toupper(unlist(sequence_data))
  sequence <- gsub("[^ATGC]", "", as.character(sequence))
  sequence_str <- paste(sequence, collapse = "")
  
  ## Find all dinucleotide repeats
  find_repeats <- function(seq_str, repeat_length) {
    motifs <- sapply(1:(nchar(seq_str) - repeat_length + 1), function(i) {
      substr(seq_str, i, i + repeat_length - 1)
    })
    repeats <- table(motifs)
    repeats[repeats > 1]  # Only keep motifs that repeat
  }
  
  di_repeats <- find_repeats(sequence_str, 2)
  di_df <- data.frame(Motif = names(di_repeats), Frequency = as.numeric(di_repeats))
  di_df <- di_df[order(-di_df$Frequency), ]
  
  ## Calculate nucleotide frequencies
  nucleotides <- table(unlist(strsplit(sequence_str, "")))
  C_freq <- ifelse("C" %in% names(nucleotides), nucleotides["C"] / sum(nucleotides), 0)
  G_freq <- ifelse("G" %in% names(nucleotides), nucleotides["G"] / sum(nucleotides), 0)
  
  ## Get CG frequency (handle case where CG might not exist)
  rhino_CG <- ifelse("CG" %in% di_df$Motif, 
                     di_df$Frequency[di_df$Motif == "CG"], 
                     0)
  
  ## Calculate O/E ratio (avoid divide by zero)
  OE_CpG <- ifelse(C_freq > 0 && G_freq > 0,
                   (rhino_CG / sum(di_df$Frequency)) / (C_freq * G_freq),
                   NA)
  
  ## Return comprehensive results
  return(list(
    OE_ratio = OE_CpG,
    CG_count = rhino_CG,
    C_freq = C_freq,
    G_freq = G_freq,
    dinucleotide_df = di_df,
    sequence_length = nchar(sequence_str)
  ))
}

## Apply the function to all sequences in virus_sequences list
results <- lapply(virus_sequences, analyze_CpG_OE)


#################################################################
################### TA O/E among all sequences ###################
#################################################################
library(seqinr)
library(ggplot2)
setwd('/Users/flickazhang/Desktop/Biohackthon Project/')

# Modified analysis function for TA dinucleotides
analyze_TA_OE <- function(virus_sequence) {
  ## Extract and prepare the sequence
  sequence <- toupper(unlist(virus_seq))
  sequence <- gsub("[^ATGC]", "", as.character(sequence))
  sequence_str <- paste(sequence, collapse = "")
  
  ## Find all dinucleotide repeats
  find_repeats <- function(seq_str, repeat_length) {
    motifs <- sapply(1:(nchar(seq_str) - repeat_length + 1), function(i) {
      substr(seq_str, i, i + repeat_length - 1)
    })
    repeats <- table(motifs)
    repeats[repeats > 1]  # Only keep motifs that repeat
  }
  
  di_repeats <- find_repeats(sequence_str, 2)
  di_df <- data.frame(Motif = names(di_repeats), Frequency = as.numeric(di_repeats))
  di_df <- di_df[order(-di_df$Frequency), ]
  
  ## Calculate nucleotide frequencies
  nucleotides <- table(unlist(strsplit(sequence_str, "")))
  A_freq <- nucleotides["A"] / sum(nucleotides)
  T_freq <- nucleotides["T"] / sum(nucleotides)
  
  ## Get TA frequency
  rhino_TA <- ifelse("TA" %in% di_df$Motif, 
                     di_df$Frequency[di_df$Motif == "TA"], 
                     0)
  
  ## Calculate TA O/E ratio
  OE_TA <- (rhino_TA / sum(di_df$Frequency)) / (A_freq * T_freq)
  
  return(list(
    OE_ratio = OE_TA,
    TA_count = rhino_TA,
    A_freq = A_freq,
    T_freq = T_freq,
    dinucleotide_df = di_df,
    sequence_length = length(sequence)
  ))
  }

## Apply the function to all sequences
results_TA <- list()
for(num in 1:16) {
  seq_name <- paste0("virus_seq_", num)
  results_TA[[seq_name]] <- analyze_TA_OE(get(seq_name))
}

# Create summary data frame for TA
summary_TA <- data.frame(
  Sequence = names(results_TA),
  OE_ratio = sapply(results_TA, function(x) x$OE_ratio),
  TA_count = sapply(results_TA, function(x) x$TA_count),
  A_freq = sapply(results_TA, function(x) x$A_freq),
  T_freq = sapply(results_TA, function(x) x$T_freq),
  Length = sapply(results_TA, function(x) x$sequence_length)
)

# Sort by TA O/E ratio
summary_TA <- summary_TA[order(summary_TA$OE_ratio), ]

# Visualization
ggplot(summary_TA, aes(x = reorder(Sequence, OE_ratio), y = OE_ratio)) +
  geom_bar(stat = "identity", fill = "#66c2a5", color = "black", width = 0.7) +
  # Add text labels with rounded values
  geom_text(aes(label = round(OE_ratio, 3)), 
            vjust = -0.5,            # Position above bars
            color = "black",         # Text color
            size = 2.5) +            # Text size
  labs(title = "TpA O/E Ratios Across Virus Sequences",
       y = "TpA O/E Ratio",
       x = "Virus Sequence") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # Adjust y-axis limits to accommodate text labels
  ylim(0, max(summary_TA$OE_ratio) * 1.1)

#################################################################
###################CpG O/E Sliding Window Plot###################
#################################################################
library(zoo)
library(ggplot2)

# First, add year information to your data (assuming you have this)
# Create a vector of years corresponding to your 16 sequences
# Replace this with your actual years
years <- c(2008,2010,2009,2010,2017,2007,2009,2016,2007,2024,2002,2023,2010,2021,2021,2015)  # Example years from 2005-2020

# Add years to your summary data
summary_df$Year <- years

# Create time-based sliding window plot
window_size <- 3  # Number of years to include in each window
step_size <- 1    # Years to advance between windows

# Calculate rolling average of O/E ratios
summary_df <- summary_df[order(summary_df$Year), ]  # Ensure chronological order

# Create sliding window data
window_data <- data.frame(
  Year = rollapply(summary_df$Year, width = window_size, by = step_size, FUN = mean),
  OE_ratio = rollapply(summary_df$OE_ratio, width = window_size, by = step_size, FUN = mean),
  n_seqs = rollapply(summary_df$OE_ratio, width = window_size, by = step_size, FUN = length)
)

# Plot the temporal evolution
ggplot(window_data, aes(x = Year, y = OE_ratio)) +
  geom_line(color = "#66c2a5", size = 1.5) +
  geom_point(aes(size = n_seqs), color = "#66c2a5") +
  labs(title = "Temporal Evolution of CpG Suppression in Rhinovirus",
       subtitle = paste(window_size, "year sliding window"),
       x = "Year",
       y = "Average CpG O/E Ratio",
       size = "Sequences per window") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  scale_size_continuous(range = c(3, 6)) +
  scale_x_continuous(breaks = seq(min(years), max(years), by = 2)) +
  geom_smooth(method = "loess", span = 0.3, color = "#fc8d62", se = FALSE)


##########################################################################
###################Rhinovirus, West Nile, and Influenza###################
##########################################################################

