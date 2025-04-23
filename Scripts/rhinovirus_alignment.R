library(Biostrings)
library(DECIPHER)

# Load sequences from a FASTA file
fasta_file <- "path/to/your/sequences.fasta"
sequences <- Biostrings::readDNAStringSet(fasta_file)

# Specify directory containing FASTA files
fasta_dir <- "/Users/oreomilk/Desktop/Rhinovirus Strains"
fasta_files <- list.files(fasta_dir, pattern = "\\.txt$", full.names = TRUE)

# Read and merge sequences
sequences <- DNAStringSet()
for (file in fasta_files) {
  seq <- readDNAStringSet(file)
  sequences <- c(sequences, seq)
}

# Remove duplicates (optional)
sequences <- unique(sequences)

# Align sequences (DECIPHER handles gaps better for divergent sequences)
aligned_sequences <- AlignSeqs(
  sequences,
  iterations = 100,  # Increase for better accuracy
  refinements = 50,
  verbose = TRUE     # Track progress
)

# Save alignment
writeXStringSet(aligned_sequences, "combined_alignment.fasta")

