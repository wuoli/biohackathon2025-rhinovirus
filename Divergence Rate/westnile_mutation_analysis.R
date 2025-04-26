library(Biostrings)
library(DECIPHER)

#------------------------ Preprocessing and Alignment------------------------

# Step 0.1: Specify directory and file paths
setwd("/Users/oreomilk/Desktop/WestNile")
fasta_file <- "/Users/oreomilk/Desktop/WestNile/westnile_seq.txt"
meta_file <- "/Users/oreomilk/Desktop/WestNile/westnile_usa_meta.tsv"
alignment_file <- "westnile_alignment.fasta"

# Step 0.2: Read fasta file and meta data
all_genomes <- readDNAStringSet(fasta_file)
metadata <- read.table(meta_file, header=TRUE, sep="\t")
strain_ids <- metadata$SequenceID

# Step 1: Align sequences (DECIPHER handles gaps better for divergent sequences)
aligned_sequences <- AlignSeqs(
  all_genomes,
  iterations = 100,  # Increase for better accuracy
  refinements = 50,
  verbose = TRUE     # Track progress
)

# Step 1.2: Save alignment
writeXStringSet(aligned_sequences, alignment_file)


#------------ Mutation Rate Calculation----------------------------------------
# Load required libraries
library(Biostrings)
library(tidyverse)
library(lubridate)

# Step 1: Load aligned sequences and metadata
alignment <- readDNAStringSet(alignment_file)
metadata <- read.csv(meta_file, sep = "\t")

# Step 2: Filter by strain IDs
names(alignment) <- sub(" .*", "", names(alignment)) 
metadata <- metadata %>% filter(SequenceID %in% strain_ids)
alignment <- alignment[names(alignment) %in% strain_ids]

# Step 3: Clean dates
metadata$date_cleaned <- ymd(metadata$date_cleaned)
metadata$year <- year(metadata$date_cleaned)

# Step 4: Choose a reference sequence (earliest year ref)
ref_id <- metadata %>% arrange(year) %>% slice(1) %>% pull(SequenceID)
ref_seq <- alignment[[ref_id]]

# Step 5: Calculate pairwise differences from reference
dna_mat <- as.matrix(alignment)

# mutation_counts <- apply(dna_mat, 1, function(row) {
#   sum(row != dna_mat[ref_id, ] & dna_mat[ref_id, ] != as.raw(0) & row != as.raw(0))
# })

mutation_counts <- apply(dna_mat, 1, function(row) {
  sum(row != dna_mat[ref_id, , drop = FALSE] &  # Keep matrix dimensions
        dna_mat[ref_id, , drop = FALSE] != "-" & 
        row != "-")
})

# Merge mutation counts with metadata
results <- metadata %>%
  filter(SequenceID %in% names(mutation_counts)) %>%
  mutate(mutations = mutation_counts[SequenceID],
         years_since = year - min(year))

# Step 6: Fit linear model to estimate mutation rate
lm_fit <- lm(mutations ~ years_since, data = results)
slope <- coef(lm_fit)["years_since"]

# Step 7: Normalize by sequence length
sequence_length <- length(ref_seq)
mutation_rate <- slope / sequence_length

cat("Estimated mutation rate per site per year:", mutation_rate, "\n")
