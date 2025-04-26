library(Biostrings)
library(DECIPHER)

#------------------------ Part 0: Preprocessing and Alignment------------------------

# Step 0.1: Specify directory containing FASTA files
setwd("/Users/oreomilk/Desktop/Biohack25")
fasta_file <- "/Users/oreomilk/Desktop/Biohack25/Code/Strains/54_rhinovirus_strains/blast2-sequences.txt"
metadata_file <- "/Users/oreomilk/Desktop/Biohack25/Code/Strains/54_rhinovirus_strains/metadata.tsv"
alignment_file <- "/Users/oreomilk/Desktop/Biohack25/Code/Strains/54_rhinovirus_strains/alignment2.fasta"

# Step 0.2: read dna sequences and metadata
all_genomes <- readDNAStringSet(fasta_file)
metadata <- read.table(metadata_file, header=TRUE, sep="\t")

# Step 0.3: filter strains used
filtered_strains <- subset(metadata, Country == "USA" & !is.na(date_cleaned) & date_cleaned > 2015)
  strain_ids <- filtered_strains$SequenceID
  # same as:
  strain_ids <- c("KY348786.1", "KY369878.1", "MK520815.1", "OQ116581.1",
                  "ON311156.1", "ON311180.1", "MZ268691.1", "MZ460691.1",
                  "OK017949.1", "OK649386.1", "ON311177.1", "ON311176.1",
                  "ON311179.1", "ON311175.1", "MW969524.1", "MW969527.1",
                  "MW969536.1", "ON311155.1", "MZ363453.1", "PP158745.1",
                  "MZ268717.1", "MZ268706.1", "MZ322915.1", "MZ363437.1",
                  "MZ427505.1", "OL961524.1", "OK017918.1", "OL961531.1",
                  "MZ629145.1", "MZ670576.1", "OK161367.1", "OK161371.1",
                  "OK161404.1", "OK539509.1", "PP194066.1", "PP187406.1",
                  "PV178322.1", "PP314216.1", "PV178503.1", "PV178487.1",
                  "PV178494.1", "PV178470.1", "PV178551.1", "PV178585.1",
                  "PV178348.1", "PV178619.1", "PV178458.1", "PV178608.1",
                  "PV178548.1", "PV178386.1", "PV178591.1", "PV178425.1",
                  "PV178222.1", "PV178207.1")

# Step 0.4: Align sequences (DECIPHER handles gaps better for divergent sequences)
aligned_sequences <- AlignSeqs(
  all_genomes,
  iterations = 100,  # Increase for better accuracy
  refinements = 50,
  verbose = TRUE     # Track progress
)

# Step 0.5: Save alignment
writeXStringSet(aligned_sequences, alignment_file)

#------------------------ Part 1: Phylogenetic Tree------------------------
library(ape)
library(phangorn)

# For BioStrings: 
#alignment <- readDNAStringSet("alignment2.fasta", format = "fasta")
#base_names <- sapply(strsplit(names(alignment), " "), `[`, 1)
#names(alignment) <- base_names
#filtered_alignment <- alignment[names(alignment) %in% filtered_strains$SequenceID]

# Step 1.1: Filter Strains (metadata and alignment)
metadata <- filtered_strains
alignment <- read.dna(alignment_file, format = "fasta")
rownames(alignment) <- sapply(strsplit(rownames(alignment), " "), `[`, 1)
filtered_alignment <- alignment[rownames(alignment) %in% strain_ids, ]

# Step 1.2: Create Distance Matrix
alignment_phyDat <- phyDat(filtered_alignment, type = "DNA")
dist_matrix <- dist.ml(alignment_phyDat)

# Step 1.3: Build tree model
# Step 1.3.1: Build Neighbor Joining tree, distance-based method
tree_NJ <- NJ(dist_matrix)
#plot(tree_NJ, tip.color = tip_colors_vec, cex = 0.8, main = "Phylogenetic Tree Colored by Year")
#legend("bottomright", legend = unique(dates$year), 
#       col = viridis(length(unique(dates$year))), 
#       pch = 19, title = "Year", cex = 0.6)

# Step 1.3.2: Build Maximum Likelihood Tree
fit <- pml(tree_NJ, data = alignment_phyDat)
fit_opt <- optim.pml(fit, model = "GTR", rearrangement = "stochastic")
plot_MLT <- plot(fit_opt$tree, main = "Maximum Likelihood Tree")

# Step 1.4: Get Colors
metadata$year <- year(metadata$date_cleaned)
unique_years <- sort(unique(metadata$year))
num_years <- length(unique_years)

# Step 1.4-1: Color Option #1
#library(viridis)
#colors <- viridis(n = length(metadata$year), option = "D")[rank(metadata$year)]
#tip_colors <- setNames(colors, metadata$SequenceID)

# Step 1.4-2: Color Option #2
#palette_colors <- colorRampPalette(brewer.pal(7, "Paired"))(num_years)
#year_colors <- setNames(palette_colors, unique_years)
#tip_colors <- year_colors[as.character(metadata$year)]
#names(tip_colors) <- metadata$SequenceID

#Step 1.4-3: Color Option #3
# library(RColorBrewer)
# # Get unique years and assign colors
# unique_years <- sort(unique(metadata$year))
# num_years <- length(unique_years)
# palette_colors <- colorRampPalette(brewer.pal(7, "Paired"))(num_years)
# year_colors <- setNames(palette_colors, unique_years)
# # Assign tip colors
# tip_colors <- year_colors[as.character(metadata$year)]
# names(tip_colors) <- metadata$SequenceID

# Step 1.4-4: Color Option #4
# tip_df <- data.frame(
#   label = metadata$SequenceID,
#   year = as.factor(metadata$year)
# )
# palette_colors <- c("#0072B2", "#56B4E9", "#009E73", "#E69F00", "#D55E00", "#CC79A7", "#6F4E7C")
# names(palette_colors) <- sort(levels(tip_df$year))
# fit_opt$tree$edge.length <- fit_opt$tree$edge.length * 0.8  # Shorter branches

# Step 1.4-5: Color Option #5
library(RColorBrewer)
tip_df <- data.frame(
  label = metadata$SequenceID,
  year = as.factor(metadata$year)
)
green_palette <- colorRampPalette(c("black", "#66c2a5"))(7)
names(green_palette) <- sort(levels(tip_df$year))
fit_opt$tree$edge.length <- fit_opt$tree$edge.length * 0.8  # Shorter branches

# Step 1.5: Plot Tree
library(ggtree)
library(treeio)
library(ggplot2)
# Plot
p <- ggtree(fit_opt$tree, layout = "rectangular") %<+% tip_df +
  geom_tiplab(aes(color = year), size = 3.2, fontface = "bold.italic", hjust = -0.1, show.legend = FALSE) +
  scale_color_manual(values = green_palette, name = "Year") +  # Set legend title here
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

saveRDS(tree_NJ, file = "tree_NJ.rds")
saveRDS(fit_opt, file = "fit_opt.rds")

#Bootstrapping
bs <- bootstrap.pml(fit_opt, bs = 100, optNni = TRUE)
plotBS(fit_opt$tree, bs, p = 50)

#------------------------ Part 2: Mutation Rate Calculation------------------------
# Load required libraries
library(Biostrings)
library(tidyverse)
library(lubridate)

# Step 3.1: Filter Strains 
metadata$year <- year(ymd(metadata$date_cleaned))
alignment <- readDNAStringSet(alignment_file)
names(alignment) <- sapply(strsplit(names(alignment), " "), `[`, 1)
filtered_alignment <- alignment[names(alignment) %in% filtered_strains$SequenceID]
dna_mat <- as.matrix(filtered_alignment)

# Step 3.2: Get ref seq
ref_id <- metadata %>% arrange(year) %>% slice(1) %>% pull(SequenceID)
ref_seq <- alignment[[ref_id]]

# Step 3.3: Calculate pairwise differences from reference
# calc_mutations <- function(seq) {
#   sum(as.character(seq) != as.character(ref_seq) &
#         as.character(seq) != "-" &
#         as.character(ref_seq) != "-")
# }
# mutation_counts <- sapply(alignment[names(alignment) %in% metadata$SequenceID], calc_mutations)

# Count the differences between current sequence and reference sequence
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

# Step 6-1: Fit linear model to estimate mutation rate by sequence length
lm_fit <- lm(mutations ~ years_since, data = results)
slope <- coef(lm_fit)["years_since"]
sequence_length <- length(ref_seq)
mutation_rate <- slope / sequence_length

# Step 6-2: Estimate mutation rate per base for each strand compared to reference
#mutation_rate2 <- results$mutations/results$years_since/sequence_length
#mean(mutation_rate2[-c(1:3)], na.rm=TRUE)

cat("Estimated mutation rate per site per year:", mutation_rate, "\n")

write.table(results, "mutations.txt", sep="\t", row.names=FALSE, quote=FALSE )

saveRDS(lm_fit, file = "mutation_lm_fit.rds")

#------------------------ Part 3: Mutation Rate Calculation for Some Regions -----------
library(Biostrings)
library(tidyverse)
library(lubridate)

# Step 3.1: Filter Strains 
metadata$year <- year(ymd(metadata$date_cleaned))
alignment <- readDNAStringSet(alignment_file)
names(alignment) <- sapply(strsplit(names(alignment), " "), `[`, 1)
filtered_alignment <- alignment[names(alignment) %in% filtered_strains$SequenceID]
dna_mat <- as.matrix(filtered_alignment)

# Step 3.2: Get ref seq
ref_id <- metadata %>% arrange(year) %>% slice(1) %>% pull(SequenceID)
ref_seq <- alignment[[ref_id]]

# Step 3: Define genomic features/regions
features <- data.frame(
  feature = c("Start Non-coding", "Genome Polyprotein", "Capsid Protein VP4", "Capsid Protein VP2", "Capsid Protein VP3", 
              "Capsid Protein VP1", "Protease 2A", "Protein 2B", "Protein 2C", "Protein 3A", "Protein 3B", 
              "Protease 3C", "Protein 3D", "End Non-coding"),
  start = c(0, 616, 616, 817, 1600, 2305, 3127, 3553, 3850, 4828, 5053, 5119, 5668, 7051),
  end = c(615, 7050, 816, 1599, 2304, 3126, 3552, 3849, 4827, 5052, 5118, 5667, 7047, length(ref_seq))
)

# Step 4: get average mutation rate per region
## Initialize results list
region_results <- list()

## Loop through each region
for (i in seq_len(nrow(features))) {
  region <- features[i, ]
  region_name <- region$feature
  region_start <- region$start
  region_end <- region$end
  region_len <- region_end - region_start + 1
  
  # Slice DNA matrix for region
  region_mat <- dna_mat[, region_start:region_end]
  
  # Calculate mutation counts per sequence in this region
  ref_row <- region_mat[ref_id, ]
  mutation_counts <- apply(region_mat, 1, function(row) {
    sum(row != ref_row & ref_row != as.raw(0) & row != as.raw(0))
  })
  
  # Merge with metadata
  region_data <- metadata %>%
    filter(SequenceID %in% names(mutation_counts)) %>%
    mutate(
      mutations = mutation_counts[SequenceID],
      years_since = year - min(year)
    )
  
  # Fit linear model
  lm_fit <- lm(mutations ~ years_since, data = region_data)
  slope <- coef(lm_fit)["years_since"]
  mutation_rate <- slope / region_len
  
  # Store results
  region_results[[region_name]] <- tibble(
    feature = region_name,
    start = region_start,
    end = region_end,
    region_length = region_len,
    avg_mutation_rate = mutation_rate
  )
}

# Combine all region results
region_summary <- bind_rows(region_results)

# View results
print(region_summary %>% arrange(desc(avg_mutation_rate)))

# Plot 
library(ggplot2)

# Assuming you have a data frame `region_summary` with columns:
#   feature, avg_mutation_rate

feature_order = features$feature
region_summary$feature <- factor(region_summary$feature, levels = feature_order)

# by mutation rate:        x = reorder(feature, avg_mutation_rate)
ggplot(region_summary, aes(x = factor(feature, levels = rev(feature_order)), y = avg_mutation_rate)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = "Genomic Region",
    y = "Mutation Rate (subs/site/year)",
    title = "Average Mutation Rate by Viral Protein Region"
  ) +
  theme_minimal(base_size = 12)

#------------------------ Part 4: UMAP -----------
library(uwot)
library(dplyr)
library(ggplot2)

# Step 4.1-1: Convert DNAbin into character matrix
# dna_mat <- as.character(filtered_alignment)
# base_levels <- c("a","c","g","t","-", "r", "y", "m", "n", "k", "w")
# dna_factor <- apply(dna_mat, 2, function(col) factor(col, levels = base_levels))
# #Chehck NAs: dna_mat[is.na(dna_factor) ]
# one_hot_encoded <- model.matrix(~ . - 1, as.data.frame(dna_factor))
# num_mat <- matrix(match(tolower(dna_mat), base_levels), 
#                   nrow = nrow(dna_mat), 
#                   ncol = ncol(dna_mat))

# Step 4.0: Filter Strains (metadata and alignment)
metadata <- filtered_strains
alignment <- read.dna(alignment_file, format = "fasta")
rownames(alignment) <- sapply(strsplit(rownames(alignment), " "), `[`, 1)
filtered_alignment <- alignment[rownames(alignment) %in% strain_ids, ]

# Step 4.1-2: Create Distance Matrix
dist_mat <- dist.dna(filtered_alignment, model = "N", pairwise.deletion = TRUE)
#library(umap)
#umap_result <- umap(as.matrix(dist_mat))
#umap_df <- as.data.frame(umap_result)
#colnames(umap_df) <- c("UMAP1", "UMAP2")
#umap_df$Strain <- rownames(filtered_alignment) # add strain names

library(uwot)
umap_result <- uwot::umap(dist_mat,
                    n_neighbors = 10,    # Try 5 to 50
                    min_dist = 0.5,     # Try 0.001 to 0.5
                    metric = "euclidean")
umap_df <- as.data.frame(umap_result)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Strain <- rownames(filtered_alignment) # add strain names

#umap_result <- umap(as.matrix(dist_mat), config = umap_config)

# Step 4.2: Add Year
metadata$year <- year(metadata$date_cleaned)
umap_df <- umap_df %>%
  left_join(metadata, by = c("Strain" = "SequenceID"))

# Step 4.2: Plot
#palette_colors <- c("#0072B2", "#56B4E9", "#009E73", "#E69F00", "#D55E00", "#CC79A7", "#6F4E7C")
#red_palette <- colorRampPalette(c("black", "red"))(7)
green_palette <- colorRampPalette(c("black", "#66c2a5"))(7)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = as.factor(year))) +
  geom_point(size = 2) +
  theme_minimal() +
  scale_color_manual(values = green_palette) +
  ggtitle("UMAP of Rhinovirus Strains Colored by Year") +
  labs(color = "Year") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

write.table(umap_df, "umap.txt", sep="\t", row.names=FALSE, quote=FALSE )

#------------------------ Part 5: BarPLot Comparison-----
data <- data.frame(
  Virus = c("Rhinovirus", "West Nile"),
  #MutationRate = c(0.01333138, 0.000105) # values if - is counted 
  MutationRate = c(0.01218725, 0.0001079434) # _ is not counted 
)

# Create the barplot
ggplot(data, aes(x = Virus, y = MutationRate, fill = Virus)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(
    y = "Sequence Divergence Rate",
    x = "Virus",
    title = "Sequence Divergence Rate of Viruses"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, color="black"),
    legend.position = "none"
    
  ) + scale_fill_manual(values = c("#66c2a5", "#8da0cb"))

#-------UMAP (Discard)-----

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
  n_neighbors = 10,
  min_dist = 0.5
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
