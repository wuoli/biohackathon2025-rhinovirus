library(seqinr)
library(ggplot2)
library(Biostrings)
library(dplyr)
# Set your working directory to where you want to save the file (adjust the path as needed)
setwd("/Users/priya/Downloads/Biohacks")  

virus_seq <- readDNAStringSet("rhino_1.txt")[[1]]  # The [[1]] loads the first sequence

# Convert to character vector of single bases
seq_chars <- unlist(strsplit(as.character(virus_seq), split = ""))

# Make sure sequence length is divisible by 3 (trim if necessary)
trim_len <- length(seq_chars) - length(seq_chars) %% 3
seq_trimmed <- seq_chars[1:trim_len]

# Group into codons
codons <- sapply(seq(1, length(seq_trimmed), by = 3), function(i) {
  paste(seq_trimmed[i:(i+2)], collapse = "")
})

# Count codon frequencies
codon_table <- table(codons)
codon_freq <- codon_table / sum(codon_table)

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

human_codon_freq <- human_codon_usage / sum(human_codon_usage)
common_codons <- intersect(names(codon_freq), names(human_codon_freq))
codon_freq_aligned <- codon_freq[common_codons]
human_codon_freq_aligned <- human_codon_freq[common_codons]

f_hum <- max(human_codon_freq[c("TTT", "TTC")])
l_hum <- max(human_codon_freq[c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG")])
s_hum <- max(human_codon_freq[c("TCT", "TCC", "TCA", "TCG", "AGT", "AGC")])
y_hum <- max(human_codon_freq[c("TAT", "TAC")])
stop_hum <- max(human_codon_freq[c("TAA", "TAG", "TGA")])
c_hum <- max(human_codon_freq[c("TGT", "TGC")])
w_hum <- human_codon_freq["TGG"]
p_hum <- max(human_codon_freq[c("CCT", "CCC", "CCA", "CCG")])
h_hum <- max(human_codon_freq[c("CAT", "CAC")])
q_hum <- max(human_codon_freq[c("CAA", "CAG")])
r_hum <- max(human_codon_freq[c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG")])
i_hum <- max(human_codon_freq[c("ATT", "ATC", "ATA")])
m_hum <- human_codon_freq["ATG"]
t_hum <- max(human_codon_freq[c("ACT", "ACC", "ACA", "ACG")])
n_hum <- max(human_codon_freq[c("AAT", "AAC")])
k_hum <- max(human_codon_freq[c("AAA", "AAG")])
v_hum <- max(human_codon_freq[c("GTT", "GTC", "GTA", "GTG")])
a_hum <- max(human_codon_freq[c("GCT", "GCC", "GCA", "GCG")])
d_hum <- max(human_codon_freq[c("GAT", "GAC")])
e_hum <- max(human_codon_freq[c("GAA", "GAG")])
g_hum <- max(human_codon_freq[c("GGT", "GGC", "GGA", "GGG")])

w_TTT <- codon_freq["TTT"] / f_hum
w_TTC <- codon_freq["TTC"] / f_hum

w_TTA <- codon_freq["TTA"] / l_hum
w_TTG <- codon_freq["TTG"] / l_hum
w_CTT <- codon_freq["CTT"] / l_hum
w_CTC <- codon_freq["CTC"] / l_hum
w_CTA <- codon_freq["CTA"] / l_hum
w_CTG <- codon_freq["CTG"] / l_hum

w_TCT <- codon_freq["TCT"] / s_hum
w_TCC <- codon_freq["TCC"] / s_hum
w_TCA <- codon_freq["TCA"] / s_hum
w_TCG <- codon_freq["TCG"] / s_hum
w_AGT <- codon_freq["AGT"] / s_hum
w_AGC <- codon_freq["AGC"] / s_hum

w_TAT <- codon_freq["TAT"] / y_hum
w_TAC <- codon_freq["TAC"] / y_hum

w_TAA <- codon_freq["TAA"] / stop_hum
w_TAG <- codon_freq["TAG"] / stop_hum
w_TGA <- codon_freq["TGA"] / stop_hum

w_TGT <- codon_freq["TGT"] / c_hum
w_TGC <- codon_freq["TGC"] / c_hum

w_TGG <- codon_freq["TGG"] / w_hum

w_CCT <- codon_freq["CCT"] / p_hum
w_CCC <- codon_freq["CCC"] / p_hum
w_CCA <- codon_freq["CCA"] / p_hum
w_CCG <- codon_freq["CCG"] / p_hum

w_CAT <- codon_freq["CAT"] / h_hum
w_CAC <- codon_freq["CAC"] / h_hum

w_CAA <- codon_freq["CAA"] / q_hum
w_CAG <- codon_freq["CAG"] / q_hum

w_CGT <- codon_freq["CGT"] / r_hum
w_CGC <- codon_freq["CGC"] / r_hum
w_CGA <- codon_freq["CGA"] / r_hum
w_CGG <- codon_freq["CGG"] / r_hum
w_AGA <- codon_freq["AGA"] / r_hum
w_AGG <- codon_freq["AGG"] / r_hum

w_ATT <- codon_freq["ATT"] / i_hum
w_ATC <- codon_freq["ATC"] / i_hum
w_ATA <- codon_freq["ATA"] / i_hum

w_ATG <- codon_freq["ATG"] / m_hum

w_ACT <- codon_freq["ACT"] / t_hum
w_ACC <- codon_freq["ACC"] / t_hum
w_ACA <- codon_freq["ACA"] / t_hum
w_ACG <- codon_freq["ACG"] / t_hum

w_AAT <- codon_freq["AAT"] / n_hum
w_AAC <- codon_freq["AAC"] / n_hum

w_AAA <- codon_freq["AAA"] / k_hum
w_AAG <- codon_freq["AAG"] / k_hum

w_GTT <- codon_freq["GTT"] / v_hum
w_GTC <- codon_freq["GTC"] / v_hum
w_GTA <- codon_freq["GTA"] / v_hum
w_GTG <- codon_freq["GTG"] / v_hum

w_GCT <- codon_freq["GCT"] / a_hum
w_GCC <- codon_freq["GCC"] / a_hum
w_GCA <- codon_freq["GCA"] / a_hum
w_GCG <- codon_freq["GCG"] / a_hum

w_GAT <- codon_freq["GAT"] / d_hum
w_GAC <- codon_freq["GAC"] / d_hum

w_GAA <- codon_freq["GAA"] / e_hum
w_GAG <- codon_freq["GAG"] / e_hum

w_GGT <- codon_freq["GGT"] / g_hum
w_GGC <- codon_freq["GGC"] / g_hum
w_GGA <- codon_freq["GGA"] / g_hum
w_GGG <- codon_freq["GGG"] / g_hum

cai <- (w_TTT * w_TTC *w_TTA * w_TTG * w_CTT * w_CTC * w_CTA * w_CTG  *
  w_TAT * w_TAC * w_TCT * w_TCC * w_TCA * w_TCG * w_AGT * w_AGC *
  w_TAA * w_TAG * w_TGA * w_TGT * w_TGC * w_TGG * w_CCT *w_CCC *
  w_CCA * w_CCG * w_CAT *w_CAC * w_CAA * w_CAG * w_CGT * w_CGC *
  w_CGA * w_CGG * w_AGA * w_AGG * w_ATT * w_ATC * w_ATA * w_ATG *
  w_AAA * w_AAG * w_AAC * w_AAT * w_ACT * w_ACC * w_ACA * w_ACG *
  w_GTT * w_GTC * w_GTA * w_GTG * w_GCT * w_GCC * w_GCA * w_GCG * 
  w_GAT* w_GAC * w_GAA *w_GAG *w_GGT * w_GGC * w_GGA *w_GGG)^(1/64)

w_values <- c( w_TTT, w_TTC, w_TTA, w_TTG, w_CTT, w_CTC, w_CTA, w_CTG, 
               w_TCT, w_TCC, w_TCA, w_TCG, w_AGT, w_AGC, w_TAT, w_TAC,
               w_TAA, w_TAG, w_TGA, w_TGT, w_TGC, w_TGG, w_CCT, w_CCC, 
               w_CCA, w_CCG,w_CAT, w_CAC, w_CAA, w_CAG, w_CGT, w_CGC, 
               w_CGA, w_CGG, w_AGA, w_AGG, w_ATT, w_ATC, w_ATA, w_ATG,
              w_ACT, w_ACC, w_ACA, w_ACG, w_AAT, w_AAC, w_AAA, w_AAG,
              w_GTT, w_GTC, w_GTA, w_GTG, w_GCT, w_GCC, w_GCA, w_GCG,
             w_GAT, w_GAC, w_GAA, w_GAG, w_GGT, w_GGC, w_GGA, w_GGG
)

# Find indices of NA or 0
which(is.na(w_values) | w_values == 0)

valid_w_values <- w_values[!is.na(w_values) & w_values > 0]


cai <- prod(valid_w_values)^(1 / length(valid_w_values))

log_w_values <- log(valid_w_values)

log_sd <- sd(log_w_values)

sem_cai <- log_sd / sqrt(length(valid_w_values))

print(paste("CAI:", round(cai, 6), "±", round(sem_cai, 6)))

