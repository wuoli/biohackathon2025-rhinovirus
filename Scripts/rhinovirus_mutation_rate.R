library(rBEAST)
library(Biostrings) 
library(ape)
library(lubridate)
library(beautier)
library(beastier)

# Convert aligned sequences to NEXUS format
aligned_sequences <- read.fasta("combined_alignment.fasta")
write.nexus.data(aligned_sequences, file = "rhinovirus.nex", format = "dna")

# Sequence names
sequence_names <- c(
  "NC_009996.1", "JN205461.1", "EF582387.1", "JF907574.1",
  "MZ667421.2", "JX074056.1", "KJ675507.1", "KF734978.1",
  "MK989757.1", "MK989751.1", "KY348786.1", "MF806525.1",
  "MZ268694.1", "OK161404.1", "PV178494.1", "PV178425.1"
)

# File Type1: Formatted Dates 
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

# File Type2: Decimal Dates 
dates_decimal <- decimal_date(as.Date(dates$date))
dates_decimal <- data.frame(
  sequence = sequence_names,
  date = dates_decimal
)
write.table(
  dates_decimal,
  file = "decimal_dates.txt",
  sep = "\t",
  col.names = FALSE,  # Critical: no headers
  row.names = FALSE,
  quote = FALSE
)

# Generate BEAST XML file
beast_xml <- create_beast2_input(
  input_filename = "combined_alignment.fasta",
  tipdates_filename = "decimal_dates.txt", 
  clock_model = create_strict_clock_model(
    clock_rate_param = create_clock_rate_param(value = 1e-6)
  ),
  tree_prior = create_yule_tree_prior(),
  mcmc = create_mcmc(chain_length = 1e6, store_every = 1000)
)

# Save XML
writeLines(beast_xml, "beast_config.xml")
run_beast2("beast_config.xml")

# Step 2: Calculate decimal

log_data <- parse_beast_log("beast_config.log")
mutation_rate <- mean(log_data$clock.rate)  # Mean rate (subs/site/year)
hpd_interval <- calc_hpd(log_data$clock.rate)  # 95% confidence interval


plot_densitree(
  output$alignment_trees,
  alpha = 0.01,
  consensus = as.character(c(1:4)),
  cex = 2.0,
  scaleX = TRUE,
  scale.bar = FALSE
)


library(babette)
fasta_filename <- "combined_alignment.fasta"  
stopifnot(file.exists(fasta_filename))
output <- bbt_run(fasta_filename)

plot_densitree(
  output$alignment_trees,
  alpha = 0.01,
  consensus = as.character(c(1:4)),
  cex = 2.0,
  scaleX = TRUE,
  scale.bar = FALSE
)