library(ape)
library(phangorn)
library(lubridate)  # For date handling
library(RColorBrewer)

# Step 1: Get dates 
# Sequence names
sequence_names <- c(
  "NC_009996.1", "JN205461.1", "EF582387.1", "JF907574.1",
  "MZ667421.2", "JX074056.1", "KJ675507.1", "KF734978.1",
  "MK989757.1", "MK989751.1", "KY348786.1", "MF806525.1",
  "MZ268694.1", "OK161404.1", "PV178494.1", "PV178425.1"
)

# Formatted Dates 
formatted_dates <- c(
  "2007-06-14", "2002-01-01", "2007-04-27", "2008-09-27",
  "2009-08-24", "2009-01-01", "2010-11-07", "2010-10-20",
  "2015-07-25", "2015-11-09", "2016-01-01", "2017-01-31",
  "2021-01-01", "2021-06-15", "2023-11-24", "2024-10-22"
)
dates <- data.frame(
  sequence = sequence_names,
  date = formatted_dates
)
write.table(
  dates,
  file = "formatted_dates.txt",
  sep = "\t",
  col.names = FALSE,  # Critical: no headers
  row.names = FALSE,
  quote = FALSE
)

# Step 2: read alignment data
alignment <- read.dna("combined_alignment.fasta", format = "fasta")
alignment_phyDat <- phyDat(alignment, type = "DNA")
dist_matrix <- dist.ml(alignment_phyDat)

# Step 3: Build tree model
# Step 3, Option 1 Build a tree using a distance-based method (e.g., Neighbor Joining)
tree_NJ <- NJ(dist_matrix)
plot_NJ <- plot(tree_NJ, main = "Neighbor-Joining Tree")

# Step 3, Option 2 Using maximum likelihood
fit <- pml(tree_NJ, data = alignment_phyDat)
fit_opt <- optim.pml(fit, model = "GTR", rearrangement = "stochastic")
plot_MLT <- plot(fit_opt$tree, main = "Maximum Likelihood Tree")

# Step 5
## Step 5.1 Extract year and map to a color
dates$year <- year(dates$date)
## Step 5.2 Normalize Range to Color Scale
year_range <- range(dates$year)
normalized_years <- (dates$year - year_range[1]) / diff(year_range)

library(viridis)
colors <- viridis(n = length(dates$year))[rank(dates$year)]
# Step 5.3: Extract and clean up tip labels in the tree
tree_NJ$tip.label <- sapply(strsplit(tree_NJ$tip.label, " "), `[`, 1)
tip_colors <- setNames(colors, dates$sequence)
tip_colors_vec <- tip_colors[tree_NJ$tip.label]

# Step 6: Plot
# Plot1: Neighbor Joining Tree
plot(tree_NJ, tip.color = tip_colors_vec, cex = 0.8, main = "Phylogenetic Tree Colored by Year")
legend("bottomright", legend = unique(dates$year), 
       col = viridis(length(unique(dates$year))), 
       pch = 19, title = "Year", cex = 0.6)

# Plot2: Maximum Likelihood Tree
plot(fit_opt$tree, tip.color = tip_colors_vec, cex = 0.8, main = "Phylogenetic Tree\nMaximum Likelihood Tree")
legend("topleft", legend = unique(dates$year), 
       col = viridis(length(unique(dates$year))), 
       pch = 19, title = "Year", cex = 0.6)

#----------- Mutation Rate Calculation-----------

### Step 1: Calculate Divergence
divergence <- dist.nodes(fit_opt$tree)[1:Ntip(fit_opt$tree), Ntip(fit_opt$tree) + 1] # use ape package

### Step2: Match Dates
tip_labels <- sapply(strsplit(fit_opt$tree$tip.label, " "), `[`, 1)
sampling_years <- dates$year[match(tip_labels, dates$sequence)]

dates$year <- year(dates$date)
dates$normalized_years <- normalized_years
sampling_years <- dates$normalized_years[match(tip_labels, dates$sequence)]

### Step 3: Root-to-Tip Regression
# Linear model for mutation rate estimation
regression_model <- lm(divergence ~ sampling_years)
mutation_rate <- coef(regression_model)["sampling_years"]  

### Step 4: Visualize Temporal Signal
plot(sampling_years, divergence, 
     xlab = "Sampling Year", ylab = "Root-to-Tip Divergence (subs/site)",
     main = paste("Mutation Rate:", signif(mutation_rate, 3), "subs/site/year"),
     pch = 19, col = tip_colors_vec)
abline(regression_model, col = "red")
cat("Estimated mutation rate:", mutation_rate, "substitutions/site/year\n")
cat("R-squared:", summary(regression_model)$r.squared, "\n")

### Optional: Rate Conversion
sequence_length <- ncol(as.character(alignment_phyDat))  # Get actual length
cat("Whole-genome rate:", mutation_rate * sequence_length, "subs/genome/year\n")

#-----------Find Genes Under Positive Selection-----------


