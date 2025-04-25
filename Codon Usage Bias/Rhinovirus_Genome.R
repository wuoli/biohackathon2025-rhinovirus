library(Biostrings)
library(ggplot2)
library(dplyr)
# Define the file names (assuming they're all in the current working directory)
setwd("/Users/priya/Downloads/Biohacks")  

source_path <- "/Users/priya/Downloads/Biohacks"  

file_names <- paste0("rhino_", 1:16, ".txt")

# Initialize an empty list to store codon tables
codon_tables <- list()

# Loop through each file
for (file in file_names) {
  # Read the sequence from the current file
  all_seqs <- readDNAStringSet(file)
  
  # Process the first (and only) sequence in the file (if it's just one sequence per file)
  seq_chars <- unlist(strsplit(as.character(all_seqs[[1]]), split = ""))
  
  # Trim to make length divisible by 3
  trim_len <- length(seq_chars) - length(seq_chars) %% 3
  seq_trimmed <- seq_chars[1:trim_len]
  
  # Group into codons
  codons <- sapply(seq(1, length(seq_trimmed), by = 3), function(i) {
    paste(seq_trimmed[i:(i+2)], collapse = "")
  })
  
  # Count codon frequencies
  codon_table <- table(codons)
  codon_freq <- codon_table / sum(codon_table)
  
  # Store as named vector
  codon_tables[[file]] <- codon_freq
}

# Combine all frequency tables into a data frame

all_codons <- unique(unlist(lapply(codon_tables, names)))
codon_matrix <- sapply(codon_tables, function(freq) {
  # Fill in missing codons with 0
  freq_full <- setNames(rep(0, length(all_codons)), all_codons)
  freq_full[names(freq)] <- freq
  return(freq_full)
})

codon_means <- rowMeans(codon_matrix)
codon_sds <- apply(codon_matrix, 1, sd)
codon_means_normalized <- codon_means / sum(codon_means)
codon_sds_normalized <- codon_sds / sum(codon_means)

# Create data frame for plotting
codon_df <- data.frame(
  Codon = names(codon_means_normalized),
  Frequency = codon_means_normalized,
  SD = codon_sds_normalized
)
codon_df <- codon_df[codon_df$Codon != "NTC", ]

# Recalculate the total sum of frequencies after exclusion
total_sum <- sum(codon_df$Frequency)
codon_df$SD <- codon_df$SD / total_sum

# Normalize the frequencies to make sure they sum to 1
codon_df$Frequency <- codon_df$Frequency / total_sum
# Plot
ggplot(codon_df, aes(x = reorder(Codon, -Frequency), y = Frequency)) +
  geom_bar(stat = "identity", fill = "#66c2a5") +
  geom_errorbar(aes(ymin = Frequency - SD, ymax = Frequency + SD), width = 0.4) +
  coord_flip() + 
  theme_minimal() +
  labs(title = "Average Codon Usage Frequency in 16 Rhinovirus Strains",
       x = "Codon", y = "Normalized Frequency") +
  theme(axis.text.x = element_text(size = 10.5, angle = 90, vjust = 0.5, hjust = 1),  # Frequency axis (now horizontal)
        axis.text.y = element_text())


human_codon_usage <- c(
  # Amino acids and their codons
  "TTT" = 0.0169, "TTC" = 0.0204, "TTA" = 0.0072, "TTG" = 0.0126, 
  "TAT" = 0.0120, "TAC" = 0.0156, "TAA" = 0.0007, "TAG" = 0.0005,
  "CTT" = 0.0128, "CTC" = 0.0194, "CTA" = 0.0069, "CTG" = 0.0403,
  "CAT" = 0.0104, "CAC" = 0.0149, "CAA" = 0.0118, "CAG" = 0.0346,
  "ATT" = 0.0157, "ATC" = 0.0214, "ATA" = 0.0071, "ATG" = 0.0223,
  "AAT" = 0.0167, "AAC" = 0.0195, "AAA" = 0.0240, "AAG" = 0.0329, 
  "GTT" = 0.0109, "GTC" = 0.0146, "GTA" = 0.0070, "GTG" = 0.0289,
  "GAT" = 0.0223, "GAC" = 0.0260, "GAA" = 0.0290, "GAG" = 0.0408,
  "TCT" = 0.0146, "TCC" = 0.0174, "TCA" = 0.0117, "TCG" = 0.0045,
  "TGT" = 0.0099, "TGC" = 0.0122, "TGA" = 0.0013, "TGG" = 0.0128,
  "CCT" = 0.0173, "CCC" = 0.0200, "CCA" = 0.0167, "CCG" = 0.0070,
  "CGT" = 0.0047, "CGC" = 0.0109, "CGA" = 0.0063, "CGG" = 0.0119,
  "ACT" = 0.0128, "ACC" = 0.0192, "ACA" = 0.0148, "ACG" = 0.0062,
  "AGT" = 0.0119, "AGC" = 0.0194, "AGA" = 0.0115, "AGG" = 0.0114,
  "GCT" = 0.0186, "GCC" = 0.0285, "GCA" = 0.0160, "GCG" = 0.0076,
  "GGT" = 0.0108, "GGC" = 0.0228, "GGA" = 0.0163, "GGG" = 0.0164
) 

human_df <- data.frame(
  Codon = names(human_codon_usage),
  Frequency = as.numeric(human_codon_usage),
  Species = "Human"
)


rhino_codon_usage <- codon_means_normalized / sum(codon_means_normalized)
rhino_df <- data.frame(
  Codon = names(rhino_codon_usage),
  Frequency = as.numeric(rhino_codon_usage),
  Species = "Rhinovirus"
)


# Step 2: Combine both dataframes
combined_df <- rbind(human_df, rhino_df)

# Ensure both codon usage data are aligned (i.e., only keep codons that appear in both)
codon_frequencies <- codon_freq[names(codon_freq) %in% names(human_codon_usage)]
human_codon_usage <- human_codon_usage[names(codon_frequencies)]  # Align the human codon usage table

# Visualize differences using a bar chart
df <- data.frame(
  Codon = rep(names(codon_frequencies), 2),
  Frequency = c(codon_frequencies, human_codon_usage),
  Type = rep(c("Rhinovirus", "Human"), each = length(codon_frequencies))
)

common_codons <- intersect(names(human_codon_usage), names(codon_freq))

# For rhinovirus and human, calculate the standard deviation of their frequencies
codon_sds_normalized <- codon_sds / sum(codon_means)  # Already calculated earlier
rhino_sds <- codon_sds_normalized[common_codons]

# Calculate the difference (virus - human)
diff_df <- data.frame(
  Codon = common_codons,
  Rhino = as.numeric(rhino_codon_usage[common_codons]),
  Human = as.numeric(human_codon_usage[common_codons]),
  Rhino_SD = as.numeric(rhino_sds)
)
diff_df$Difference <- diff_df$Rhino - diff_df$Human
diff_df$Lower <- diff_df$Difference - diff_df$Rhino_SD
diff_df$Upper <- diff_df$Difference + diff_df$Rhino_SD

print(diff_df)

ggplot(diff_df, aes(x = reorder(Codon, -Difference), y = Difference, fill = Difference > 0)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.4) +
  scale_fill_manual(values = c("TRUE" = "#8dc0cd", "FALSE" = "#fc8d62"),
                    labels = c("Underused", "Overused"),
                    name = "Rhinovirus Usage") +
  theme_minimal() +
  labs(title = "Codon Usage Difference: Rhinovirus vs. Human",
       x = "Codon", y = "Difference in Frequency (Virus - Human)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12), axis.text.y = element_text(size = 12))+
  scale_x_discrete(expand = expansion(add = c(0.6, 0.6)))  # Adjust spacing)



