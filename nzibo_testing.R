coldata <- readRDS("~/OneDrive - University of Otago/GenomicsAotearoa/Conferences-meetings-engagements/NZIBO/NZIBO-workshop-april2026/data/coldata.rds")
counts <- readRDS("~/OneDrive - University of Otago/GenomicsAotearoa/Conferences-meetings-engagements/NZIBO/NZIBO-workshop-april2026/data/countdata.rds")
resultsPadjLogFC <- readRDS("~/OneDrive - University of Otago/GenomicsAotearoa/Conferences-meetings-engagements/NZIBO/NZIBO-workshop-april2026/data/resultsPadjLogFC.rds")

countdata <- read.csv("https://raw.githubusercontent.com/GenomicsAotearoa/NZIBO-bioinformatics-workshop-2026/refs/heads/main/data/count.data.csv", row.names = "GeneID")
head(countdata) 





library(DESeq2)

coldata$dex <-  relevel(factor(coldata$dex), ref = "untrt")
coldata$cell <- factor(coldata$cell)
str(coldata)


dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData   = coldata,
  design    = ~ cell + dex
)

dds <- DESeq(dds)

results <- DESeq2::results(dds)

results <- na.omit(results)


resultsNames(dds)



library(ggplot2)
library(ggiraph)
library(plotly)


# Define significance thresholds
pval_threshold <- 0.05
fc_threshold <- 1  # log2 fold change threshold

# Create a data frame for plotting
plot_data <- data.frame(
  logFC = results$log2FoldChange,
  negLogPval = -log10(results$pvalue),
  adj.P.Val = results$padj,
  ID = rownames(results)  # Assuming row names are gene IDs
)

# Add a column to categorize genes
plot_data$category <- ifelse(plot_data$adj.P.Val <= pval_threshold,
                             ifelse(plot_data$logFC >= fc_threshold, "Upregulated",
                                    ifelse(plot_data$logFC <= -fc_threshold, "Downregulated", "Passes P-value cut off")),
                             "Not Significant")


# Create the ggplot object
p <- ggplot(plot_data, aes(x = logFC, y = negLogPval, color = category, text = ID)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey20", "Passes P-value cut off" = "grey")) +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed") +
  labs(
    title = "Interactive Volcano Plot of Differential Gene Expression",
    subtitle = paste("Thresholds: |log2FC| >", fc_threshold, "and adjusted p-value <", pval_threshold),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

# Convert ggplot to an interactive plotly object
interactive_plot <- ggplotly(p, tooltip = c("text", "x", "y", "color"))

# Customize hover text
interactive_plot <- interactive_plot %>% 
  layout(hoverlabel = list(bgcolor = "white"),
         hovermode = "closest")

# Display the interactive plot
interactive_plot

# display flat plot 
p 







resultsPadjLogFC_df <- as.data.frame(resultsPadjLogFC)
resultsUpReg <- resultsPadjLogFC_df[resultsPadjLogFC_df$log2FoldChange>=0,]
resultsDownReg <- resultsPadjLogFC_df[resultsPadjLogFC_df$log2FoldChange<=0,]

#check split work
nrow(resultsUpReg) + nrow(resultsDownReg) == nrow(resultsPadjLogFC_df)
#Upreg = treated
#Downreg = untreated

head(rownames(resultsUpReg))
#"ENSG00000179593" 
#counts[counts$GeneID == "ENSG00000179593",] #if gene ID own column
counts["ENSG00000179593",] #if row.names = T

# yes confirm - up reg = treated. 


resultsUpReg <- resultsUpReg |> arrange(desc(log2FoldChange))
resultsDownReg <- resultsDownReg |> arrange(log2FoldChange)
nrow(resultsUpReg)#528
nrow(resultsDownReg)#483

# top 3 upregulated genes in treatment

head(resultsUpReg)
head(rownames(resultsUpReg), 3)

"ENSG00000179593" "ENSG00000109906" "ENSG00000250978" 

library(biomaRt)

mart <- useMart("ensembl", 
                dataset = "hsapiens_gene_ensembl",
                host = "https://useast.ensembl.org")   # US East mirror

# other mirrors to try:
# "https://uswest.ensembl.org"
# "https://asia.ensembl.org"

genes3UP <- c("ENSG00000179593" ,"ENSG00000109906" ,"ENSG00000250978")

info3UP <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "hgnc_symbol",
    "description",
    "chromosome_name",
    "go_id",
    "name_1006",
    "gene_biotype"
  ),
  filters = "ensembl_gene_id",
  values = genes3UP,
  mart = mart
)


# top 5 downregulated genes in treatment

head(resultsDownReg)
head(rownames(resultsDownReg), 5)

"ENSG00000128285" "ENSG00000267339" "ENSG00000019186" "ENSG00000183454"
"ENSG00000146006"

library(biomaRt)

#mart <- useMart("ensembl", 
                dataset = "hsapiens_gene_ensembl",
                host = "https://useast.ensembl.org")   # US East mirror

# other mirrors to try:
# "https://uswest.ensembl.org"
# "https://asia.ensembl.org"

genes5down <- c("ENSG00000128285", "ENSG00000267339", "ENSG00000019186" ,"ENSG00000183454",
           "ENSG00000146006")

infodown <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "hgnc_symbol",
    "description",
    "chromosome_name",
    "go_id",
    "name_1006",
    "gene_biotype"
  ),
  filters = "ensembl_gene_id",
  values = genes5down,
  mart = mart
)



