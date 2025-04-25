library(enrichR)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggnewscale)
library(dplyr)
setwd('/Users/flickazhang/Desktop/Biohackthon_Project/')
dbs <- listEnrichrDbs()
websiteLive <- TRUE
if (websiteLive) {
  dbs <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", 
           "GO_Molecular_Function_2023", "KEGG_2021_Human")
}
# Set up databases for your project
dbs.use <- c("Biocarta_2016",           # Biocarta pathways
             "Reactome_2016",           # Reactome pathways
             "GO_Biological_Process_2018",  # GO Biological Process
             "GO_Molecular_Function_2018",  # GO Molecular Function
             "GO_Cellular_Component_2018") # GO Cellular Component

rhino_gene_list<-read.csv("/Users/flickazhang/Desktop/Biohackthon_Project/Rhinovirus.csv")
rhino_sig_genes <- rhino_gene_list[rhino_gene_list$padj < 0.05,]
rhino_sig_genes<-rhino_sig_genes$Gene
wn_gene_list<-read.csv("/Users/flickazhang/Desktop/Biohackthon_Project/WestNile.csv")
wn_sig_genes <- wn_gene_list[wn_gene_list$padj < 0.05,]
wn_sig_genes<-wn_sig_genes$Gene

deg.list <- list(rhino_sig_genes, wn_sig_genes)
names(deg.list) <- c("Rhinovirus", "West Nile")

# Function to perform enrichment analysis and generate dot plots
perform_enrichment <- function(ontology) {
  # Perform enrichment analysis
  enrich.comp <- compareCluster(
    geneClusters = deg.list,
    fun = "enrichGO",
    ont = ontology,  # BP, CC, or MF
    OrgDb = "org.Hs.eg.db",
    keyType = "SYMBOL",
    pvalueCutoff = 0.05
  )
  p <- enrichplot::dotplot(enrich.comp, showCategory = 7) +
    theme(
      axis.text.x = element_text(size = 22, angle = 35, hjust = 1),  # Larger x-axis labels
      axis.text.y = element_text(size = 19),  # Larger y-axis labels
      axis.title.x = element_text(size = 24, face = "bold"),  # Larger x-axis title
      axis.title.y = element_text(size = 22, face = "bold"),  # Larger y-axis title
      legend.text = element_text(size = 22),  # Larger legend text
      legend.title = element_text(size = 24, face = "bold"),  # Larger legend title
      panel.spacing = unit(0.001, "lines")  # Reduce space between columns
    )
  
  # Save the plot
  ggsave(
    paste0("./enrichment_pathway_analysis_wn_rhino", ontology, ".tiff"),
    plot = p,
    height = 12,
    width = 8.3,
    units = "in"
  )
}

# Perform enrichment analysis for BP, CC, and MF
perform_enrichment("BP")  # Biological Process
perform_enrichment("CC")  # Cellular Component
perform_enrichment("MF")  # Molecular Function

#Biological Process
enrich.comp <- compareCluster(geneClusters = deg.list, fun = "enrichGO", ont = "BP", 
                              OrgDb="org.Hs.eg.db", keyType = "SYMBOL", pvalueCutoff=0.05)
enrich.comp <- compareCluster(geneClusters = gene_list, fun = "enrichGO", ont = "CC", 
                              OrgDb="org.Hs.eg.db", keyType = "SYMBOL", pvalueCutoff=0.05)
enrich.comp <- compareCluster(geneClusters = gene_list, fun = "enrichGO", ont = "MF", 
                              OrgDb="org.Hs.eg.db", keyType = "SYMBOL", pvalueCutoff=0.05)

# Define the list of keywords
keywords <- c(
  "response to virus",
  "defense response to virus",
  "negative regulation of viralgenome replication",
  "negative regulation of viral",
  "regulation of viral genomereplication",
  "regulation of viral life cycle",
  "viral genome replication",
  "antiviral innate immune response",
  "regulation of innate immune response",
  "innate immune response-activating signaling pathway"
)

# Access the enrichment result data
df <- enrich.comp@compareClusterResult

# Filter rows that match any keyword (case-insensitive, partial match)
virus_related_df <- df[Reduce(`|`, lapply(keywords, function(k) grepl(k, df$Description, ignore.case = TRUE))), ]
enrich_df <- virus_related_df[, c(1, 3, ncol(df) - 1)]


# Optionally save to CSV
write.csv(enrich_df, "enrichmen_genes.csv", row.names = FALSE)
#remove redundant pathways

p <- enrichplot::dotplot(enrich.comp, showCategory=7) + 
  theme(axis.text.x = element_text(size = 22, angle = 35, hjust = 1),  # Larger x-axis labels
        axis.text.y = element_text(size = 19),  # Larger y-axis labels
        axis.title.x = element_text(size = 24, face = "bold"),  # Larger x-axis title
        axis.title.y = element_text(size = 22, face = "bold"),  # Larger y-axis title
        legend.text = element_text(size = 22),  # Larger legend text
        legend.title = element_text(size = 24, face = "bold"),  # Larger legend title
        panel.spacing = unit(0.001, "lines"))   # Reduce space between columns

# Save the updated plot
ggsave("./enrichment_pathway_analysis_BP.tiff", height = 12, width = 8.3, units = "in")

###########KEGG Pathway Analysis################
# Convert gene symbols to Entrez IDs (required for enrichKEGG)
GCI.entrez <- bitr(GCI.genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
LB.entrez <- bitr(LB.genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
PFF.entrez <- bitr(PFF.genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")

# Perform KEGG pathway enrichment for each group
kegg.GCI <- enrichKEGG(gene = GCI.entrez$ENTREZID, organism = "mmu", pvalueCutoff = 0.05)
kegg.LB <- enrichKEGG(gene = LB.entrez$ENTREZID, organism = "mmu", pvalueCutoff = 0.05)
kegg.PFF <- enrichKEGG(gene = PFF.entrez$ENTREZID, organism = "mmu", pvalueCutoff = 0.05)

# Extract Top 10 significant pathways for each group
top10.GCI <- kegg.GCI@result %>% arrange(pvalue) %>% head(10)
top10.LB <- kegg.LB@result %>% arrange(pvalue) %>% head(10)
top10.PFF <- kegg.PFF@result %>% arrange(pvalue) %>% head(10)

write.csv(top10.GCI, "./Top10_KEGG_GCI.csv", row.names = FALSE)
write.csv(top10.LB, "./Top10_KEGG_LB.csv", row.names = FALSE)
write.csv(top10.PFF, "./Top10_KEGG_PFF.csv", row.names = FALSE)
###########################

plot <- ggplot(enrich.comp, aes(x = Cluster, y = Description, color = Cluster, size = -log10(p.adjust))) +
  geom_point() +
  scale_color_manual(values = cluster_colors) +
  scale_size_continuous(breaks = c(1.32, 5, 10, 20, 30), labels = c("-log10(0.05)", 5, 10, 20, 30), name = "-log10(FDR)") +
  labs(x = "Cluster", y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 33, hjust = 1),
        axis.title.y = element_text(size = 14),
        plot.background = element_rect(color = "white", fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.25)) +
  guides(color = FALSE)
ggsave("./enrichment_2.tiff", height = 6, width = 6, units = "in")

Collapse









