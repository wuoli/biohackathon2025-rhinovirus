
#------------------------ Part 0: Preprocessing and Alignment------------------------
# Step 0.1: Specify directory containing FASTA files
setwd("/Users/oreomilk/Desktop/Rhinovirus Strains")
fasta_file <- "/Users/oreomilk/Desktop/Rhinovirus Strains/all/sequences.fasta"
metadata_file <- "/Users/oreomilk/Desktop/Rhinovirus Strains/all/sequences.tsv"
alignment_file <- "/Users/oreomilk/Desktop/Rhinovirus Strains/all/alignment.fasta"

# Step 0.2: read dna sequences and metadata
all_genomes <- readDNAStringSet(fasta_file)
base_names <- sapply(strsplit(names(all_genomes), " "), `[`, 1)
metadata <- read.table(metadata_file, header=TRUE, sep="\t", fill=NA)
names(all_genomes) <- base_names

# Step 0.3: filter strains used
filtered_strains <- subset(metadata, Country == "USA" & !is.na(Collection_Date) & Collection_Date != "")
strain_ids <- filtered_strains$SequenceID

usa_genomes <- all_genomes[names(all_genomes) %in% filtered_strains$Accession]

# Step 0.4: Align sequences (DECIPHER handles gaps better for divergent sequences)
aligned_sequences <- AlignSeqs(
  usa_genomes,
  iterations = 100,  # Increase for better accuracy
  refinements = 50,
  verbose = TRUE     # Track progress
)

# Step 0.5: Save alignment
writeXStringSet(aligned_sequences, alignment_file)

#------------------------ Part 1: Phylogenetic Tree------------------------
library(ape)
library(phangorn)

filtered_alignment <- read.dna(alignment_file, format = "fasta")

# Step 1.2: Create Distance Matrix
alignment_phyDat <- phyDat(filtered_alignment, type = "DNA")
dist_matrix <- dist.ml(alignment_phyDat)
tree_NJ <- NJ(dist_matrix)
fit <- pml(tree_NJ, data = alignment_phyDat)
fit_opt <- optim.pml(fit, model = "GTR", rearrangement = "stochastic")
plot_MLT <- plot(fit_opt$tree, main = "Maximum Likelihood Tree")

# Step 1.4: Get Colors
metadata$year <- year(metadata$date_cleaned)
unique_years <- sort(unique(metadata$year))
num_years <- length(unique_years)
library(RColorBrewer)
tip_df <- data.frame(
  label = metadata$SequenceID,
  year = as.factor(metadata$year)
)
palette_colors <- c("#0072B2", "#56B4E9", "#009E73", "#E69F00", "#D55E00", "#CC79A7", "#6F4E7C")
names(palette_colors) <- sort(levels(tip_df$year))
fit_opt$tree$edge.length <- fit_opt$tree$edge.length * 0.8  # Shorter branches

# Step 1.5: Plot Tree
library(ggtree)
library(treeio)
library(ggplot2)
# Plot
p <- ggtree(fit_opt$tree, layout = "rectangular") %<+% tip_df +
  geom_tiplab(aes(color = year), size = 3.2, fontface = "bold.italic", hjust = -0.1, show.legend = FALSE) +
  scale_color_manual(values = palette_colors, name = "Year") +  # Set legend title here
  theme_tree2() +
  theme(
    legend.position = c(0.95, 0.05),
    legend.justification = c(1, 0),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10, face = "bold.italic"),
    plot.title = element_text(hjust = 0.5, size = 0),
    plot.margin = margin(10, 30, 10, 10),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank()
  ) +
  coord_cartesian(clip = 'off') +
  geom_tippoint(aes(color = year), size = 2)
print(p)



#------------------------ Part 4: UMAP -----------
library(uwot)
library(dplyr)
library(ggplot2)

filtered_alignment <- read.dna(alignment_file, format = "fasta")

# Step 4.1: Create Distance Matrix
dist_mat <- dist.dna(filtered_alignment, model = "N", pairwise.deletion = TRUE)
umap_result <- uwot::umap(dist_mat,
                          n_neighbors = 50,    # Try 5 to 50
                          min_dist = 0.3,     # Try 0.001 to 0.5
                          metric = "euclidean")
umap_df <- as.data.frame(umap_result)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Strain <- rownames(filtered_alignment) # add strain names

# Step 4.2: Add Year
library(lubridate)
library(dplyr)
dates <- filtered_strains$Collection_Date
filtered_strains$Parsed_Date <- case_when(
  grepl("^\\d{4}-\\d{2}-\\d{2}$", dates) ~ ymd(dates),       # full date
  grepl("^\\d{4}-\\d{2}$", dates)        ~ ym(dates),        # year-month
  grepl("^\\d{4}$", dates)               ~ ymd(paste0(dates, "-01-01")),  # year only
  TRUE                                   ~ NA_Date_          # fallback
)

filtered_strains$year <- year(filtered_strains$Parsed_Date)
umap_df <- umap_df %>%
  left_join(filtered_strains, by = c("Strain" = "Accession"))

# Step 4.2: Plot
palette_14 <- c(
  "#0072B2", # Deep Blue
  "#56B4E9", # Sky Blue
  "#009E73", # Teal
  "#00A99D", # Cyan-green
  "#1B7F5C", # Deep Green
  "#7CAE00", # Olive Green
  "#E69F00", # Orange
  "#D55E00", # Vermilion
  "#A6761D", # Brown
  "#CC79A7", # Magenta
  "#AA4C8F", # Purple-pink
  "#6F4E7C", # Deep Purple
  "#5B1A18", # Burgundy
  "#F564E3"  # Light Magenta
)


ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = as.factor(year))) +
  geom_point(size = 1) +
  theme_minimal() +
  scale_color_manual(values = palette_14) +
  ggtitle("UMAP of Rhinovirus Strains Colored by Year") +
  labs(color = "Year")

umap_df$year <- strain_years
filtered_strains[is.na(filtered_strains$year),]

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Year)) +
  geom_point(size = 2) +
  theme_minimal() +
  ggtitle("UMAP of Rhinovirus Strains by Year")


# Step 4.2: Create numeric coding (A=1, C=2, G=3, T=4, -=5)


# Convert DNAbin to distance matrix using Hamming distance
library(uwot)
library(ape)
dist_matrix <- dist.dna(filtered_alignment, model = "N", pairwise.deletion = TRUE)
dense_dist <- as.matrix(dist_matrix)
set.seed(123)
umap_results <- uwot::umap(
  dense_dist,
  metric = "precomputed",
  n_neighbors = 15,
  min_dist = 0.1
)
set.seed(123)
umap_results <- umap(
  as.matrix(dist_matrix),
  metric = "precomputed",
  n_neighbors = 15,
  min_dist = 0.1
)

# Step 4.3: UMAP
library(umap)
# Transpose to samples-as-rows format
numeric_data <- t(num_mat)
# Run UMAP with Hamming distance for categorical data
umap_results <- umap(numeric_data, 
                     metric = "hamming",
                     n_neighbors = 15,
                     min_dist = 0.1)

# Assuming 'years' is a vector containing collection years
umap_df <- data.frame(
  UMAP1 = umap_results$layout[,1],
  UMAP2 = umap_results$layout[,2],
  Year = years
)

library(ggplot2)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Year)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_viridis_c(name = "Collection Year") +
  labs(title = "UMAP Projection of Rhinovirus Genomic Sequences",
       subtitle = "Colored by Collection Year")
