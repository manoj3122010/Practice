library(dplyr)
library(readxl)
library(edgeR)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ReactomePA)
library(ggplot2)
library(EnhancedVolcano)

gse <- "~//shared/NGS/GSE168831_fpkm_genename.xlsx"

exp_data <- read_xlsx(gse)

selected_columns <- c("gene_id", "gene_name", "WT_1_A", "WT_1_B", "WT_1_C", "WT_2_A", "WT_2_B", "WT_2_C", 
                      "ALS_C1_A", "ALS_C1_B", "ALS_C1_C", "ALS_C2_A", "ALS_C2_B", "ALS_C2_C", 
                      "ALS_F1_A", "ALS_F1_B", "ALS_F1_C", "ALS_F2_A", "ALS_F2_B", "ALS_F2_C")

head(exp_data)

data <- exp_data[, selected_columns] 

colnames(data)

head(data)

mito_genes <- data[grep("^MT-", data$gene_name), ]
ribo_genes <- data[grep("^RPS|^RPL", data$gene_name), ]

total_expression <- colSums(data[,-c(1,2)])

total_mito_expression <- colSums(mito_genes[,-c(1,2)])
total_ribo_expression <- colSums(ribo_genes[,-c(1,2)])

percent_mito <- (total_mito_expression / total_expression) * 100
percent_ribo <- (total_ribo_expression / total_expression) * 100

percentages <- data.frame(Sample = colnames(data)[-c(1,2)],
                          Percent_Mito = percent_mito,
                          Percent_Ribo = percent_ribo)

print(percentages)

Condition <- factor(c(rep("WT", 6), rep("ALS_C", 6), rep("ALS_F", 6)))
RNA_data <- DGEList(counts=data, group=Condition)
keep <- rowSums(cpm(RNA_data) > 1) >= 2
RNA_data_filtered <- RNA_data[keep,,keep.lib.sizes=FALSE]

RNA_data_filtered <- calcNormFactors(RNA_data_filtered)
RNA_data_filtered <- estimateCommonDisp(RNA_data_filtered, verbose=TRUE)
RNA_data_filtered <- estimateTagwiseDisp(RNA_data_filtered)

colnames(RNA_data_filtered)

plotBCV(RNA_data_filtered)

et_WT_vs_ALS_C <- exactTest(RNA_data_filtered, pair=c("WT", "ALS_C"))
topTags(et_WT_vs_ALS_C)

et_WT_vs_ALS_F <- exactTest(RNA_data_filtered, pair=c("WT", "ALS_F"))
topTags(et_WT_vs_ALS_F)

et_ALS_C_vs_ALS_F <- exactTest(RNA_data_filtered, pair=c("ALS_C", "ALS_F"))
topTags(et_ALS_C_vs_ALS_F)

WTvsC9 <- et_WT_vs_ALS_C
WTvsFUS <- et_WT_vs_ALS_F
C9vsFUS <- et_ALS_C_vs_ALS_F

WTvsC9$table$gene_id <- WTvsC9$genes[,1]
WTvsFUS$table$gene_id <- WTvsFUS$genes[,1]
C9vsFUS$table$gene_id <- C9vsFUS$genes[,1]

WTvsC9$table$gene_name <- exp_data$gene_name[match(WTvsC9$table$gene_id, exp_data$gene_id)]
WTvsFUS$table$gene_name <- exp_data$gene_name[match(WTvsFUS$table$gene_id, exp_data$gene_id)]
C9vsFUS$table$gene_name <- exp_data$gene_name[match(WTvsC9$table$gene_id, exp_data$gene_id)]

WTvsC9$table$gene_description <- exp_data$gene_description[match(WTvsC9$table$gene_id, exp_data$gene_id)]
WTvsFUS$table$gene_description <- exp_data$gene_description[match(WTvsFUS$table$gene_id, exp_data$gene_id)]
C9vsFUS$table$gene_description <- exp_data$gene_description[match(WTvsC9$table$gene_id, exp_data$gene_id)]

head(WTvsC9$table)
head(WTvsFUS$table)
head(C9vsFUS$table)

write.csv(WTvsC9$table, "~//shared/NGS/WTvsC9.csv")
write.csv(WTvsFUS$table, "~//shared/NGS/WTvsFUS.csv")
write.csv(C9vsFUS$table, "~//shared/NGS/C9vsFUS.csv")

EnhancedVolcano(WTvsC9$table,
                lab = WTvsC9$table$gene_name,
                x = 'logFC',
                y = 'PValue',
                title = 'Volcano Plot: WT vs C9',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 3.5,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0)

EnhancedVolcano(WTvsFUS$table,
                lab = WTvsFUS$table$gene_name,
                x = 'logFC',
                y = 'PValue',
                title = 'Volcano Plot: WT vs FUS',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 3.5,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0)

EnhancedVolcano(C9vsFUS$table,
                lab = C9vsFUS$table$gene_name,
                x = 'logFC',
                y = 'PValue',
                title = 'Volcano Plot: C9 vs FUS',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 3.5,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0)

top_genes_WTvsC9 <- WTvsC9$table$gene_name[order(WTvsC9$table$PValue)][1:25]
top_genes_WTvsFUS <- WTvsFUS$table$gene_name[order(WTvsFUS$table$PValue)][1:25]
top_genes_C9vsFUS <- C9vsFUS$table$gene_name[order(C9vsFUS$table$PValue)][1:25]

heatmap_data_WTvsC9 <- exp_data[exp_data$gene_name %in% top_genes_WTvsC9, c('WT_1_A', 'WT_1_B', 'WT_1_C', 'WT_2_A', 'WT_2_B', 'WT_2_C', 
                                                                    'ALS_C1_A', 'ALS_C1_B', 'ALS_C1_C', 'ALS_C2_A', 'ALS_C2_B', 'ALS_C2_C')]

heatmap_data_WTvsFUS <- exp_data[exp_data$gene_name %in% top_genes_WTvsFUS, c('WT_1_A', 'WT_1_B', 'WT_1_C', 'WT_2_A', 'WT_2_B', 'WT_2_C', 
                                                                      'ALS_F1_A', 'ALS_F1_B', 'ALS_F1_C', 'ALS_F2_A', 'ALS_F2_B', 'ALS_F2_C')]

heatmap_data_C9vsFUS <- exp_data[exp_data$gene_name %in% top_genes_C9vsFUS, c('ALS_C1_A', 'ALS_C1_B', 'ALS_C1_C', 'ALS_C2_A', 'ALS_C2_B', 'ALS_C2_C', 
                                                                      'ALS_F1_A', 'ALS_F1_B', 'ALS_F1_C', 'ALS_F2_A', 'ALS_F2_B', 'ALS_F2_C')]


heatmap_data_WTvsC9 <- log2(heatmap_data_WTvsC9 + 1)
heatmap_data_WTvsFUS <- log2(heatmap_data_WTvsFUS + 1)
heatmap_data_C9vsFUS <- log2(heatmap_data_C9vsFUS + 1)

rownames(heatmap_data_WTvsC9) <- exp_data$gene_name[exp_data$gene_name %in% top_genes_WTvsC9]
rownames(heatmap_data_WTvsFUS) <- exp_data$gene_name[exp_data$gene_name %in% top_genes_WTvsFUS]
rownames(heatmap_data_C9vsFUS) <- exp_data$gene_name[exp_data$gene_name %in% top_genes_C9vsFUS]

pheatmap(heatmap_data_WTvsC9, cluster_rows=TRUE, cluster_cols=TRUE, scale="row",
         main="Heatmap: WT vs ALS_C")

pheatmap(heatmap_data_WTvsFUS, cluster_rows=TRUE, cluster_cols=TRUE, scale="row",
         main="Heatmap: WT vs ALS_F")

pheatmap(heatmap_data_C9vsFUS, cluster_rows=TRUE, cluster_cols=TRUE, scale="row",
         main="Heatmap: ALS_C vs ALS_F")

# Select top 50 genes based on P-value
top_genes_WTvsC9 <- WTvsC9$table$gene_name[order(WTvsC9$table$PValue)][1:25]
top_genes_WTvsFUS <- WTvsFUS$table$gene_name[order(WTvsFUS$table$PValue)][1:25]
top_genes_C9vsFUS <- C9vsFUS$table$gene_name[order(C9vsFUS$table$PValue)][1:25]

# Create a list of unique top genes
unique_genes <- unique(c(top_genes_WTvsC9, top_genes_WTvsFUS, top_genes_C9vsFUS))

# Create a matrix of logFC values for these genes across all comparisons
combined_logFC <- data.frame(
  WTvsC9 = WTvsC9$table$logFC[match(unique_genes, WTvsC9$table$gene_name)],
  WTvsFUS = WTvsFUS$table$logFC[match(unique_genes, WTvsFUS$table$gene_name)],
  C9vsFUS = C9vsFUS$table$logFC[match(unique_genes, C9vsFUS$table$gene_name)]
)

# Assign row names as gene names
rownames(combined_logFC) <- unique_genes

# Remove any rows with NAs (if any)
combined_logFC <- combined_logFC[complete.cases(combined_logFC), ]

# Plot the heatmap
pheatmap(
  mat = as.matrix(combined_logFC),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Heatmap of LogFC for Top Genes Across Comparisons"
)

ranked_genes_WTvsC9 <- WTvsC9$table$logFC
names(ranked_genes_WTvsC9) <- WTvsC9$table$gene_name
ranked_genes_WTvsC9 <- sort(ranked_genes_WTvsC9, decreasing=TRUE)

ranked_genes_WTvsFUS <- WTvsFUS$table$logFC
names(ranked_genes_WTvsFUS) <- WTvsFUS$table$gene_name
ranked_genes_WTvsFUS <- sort(ranked_genes_WTvsFUS, decreasing=TRUE)

ranked_genes_C9vsFUS <- C9vsFUS$table$logFC
names(ranked_genes_C9vsFUS) <- C9vsFUS$table$gene_name
ranked_genes_C9vsFUS <- sort(ranked_genes_C9vsFUS, decreasing=TRUE)

gsea_WTvsC9 <- gseGO(geneList=ranked_genes_WTvsC9, 
                     ont="BP", # Biological Process
                     keyType="SYMBOL", 
                     minGSSize=10, 
                     maxGSSize=500, 
                     pvalueCutoff=0.05, 
                     OrgDb=org.Hs.eg.db, 
                     verbose=FALSE,
                     eps = 0)

gsea_WTvsFUS <- gseGO(geneList=ranked_genes_WTvsFUS, 
                      ont="BP", 
                      keyType="SYMBOL", 
                      minGSSize=10, 
                      maxGSSize=500, 
                      pvalueCutoff=0.05, 
                      OrgDb=org.Hs.eg.db, 
                      verbose=FALSE,
                      eps = 0)

gsea_C9vsFUS <- gseGO(geneList=ranked_genes_C9vsFUS, 
                      ont="BP", 
                      keyType="SYMBOL", 
                      minGSSize=10, 
                      maxGSSize=500, 
                      pvalueCutoff=0.05, 
                      OrgDb=org.Hs.eg.db, 
                      verbose=FALSE,
                      eps = 0)

dotplot(gsea_WTvsC9, showCategory=10, title="GSEA: WT vs ALS_C")
dotplot(gsea_WTvsFUS, showCategory=10, title="GSEA: WT vs ALS_F")
dotplot(gsea_C9vsFUS, showCategory=10, title="GSEA: ALS_C vs ALS_F")

gsea_WTvsC9 <- pairwise_termsim(gsea_WTvsC9)
gsea_WTvsFUS <- pairwise_termsim(gsea_WTvsFUS)
gsea_C9vsFUS <- pairwise_termsim(gsea_C9vsFUS)

emapplot(gsea_WTvsC9)
emapplot(gsea_WTvsFUS)
emapplot(gsea_C9vsFUS)

reactome_genes_WTvsC9 <- WTvsC9$table$logFC
names(reactome_genes_WTvsC9) <- WTvsC9$table$gene_id
reactome_genes_WTvsC9 <- sort(reactome_genes_WTvsC9, decreasing=TRUE)

reactome_genes_WTvsFUS <- WTvsFUS$table$logFC
names(reactome_genes_WTvsFUS) <- WTvsFUS$table$gene_id
reactome_genes_WTvsFUS <- sort(reactome_genes_WTvsFUS, decreasing=TRUE)

reactome_genes_C9vsFUS <- C9vsFUS$table$logFC
names(reactome_genes_C9vsFUS) <- C9vsFUS$table$gene_id
reactome_genes_C9vsFUS <- sort(reactome_genes_C9vsFUS, decreasing=TRUE)

ensembl_WTvsC9 <- names(reactome_genes_WTvsC9)
ensembl_WTvsFUS <- names(reactome_genes_WTvsFUS)
ensembl_C9vsFUS <- names(reactome_genes_C9vsFUS)

head(ensembl_C9vsFUS)
head(ensembl_WTvsC9)
head(ensembl_WTvsFUS)

converted_WTvsC9 <- bitr(ensembl_WTvsC9, fromType = "ENSEMBL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)
converted_WTvsFUS <- bitr(ensembl_WTvsFUS, fromType = "ENSEMBL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)
converted_C9vsFUS <- bitr(ensembl_C9vsFUS, fromType = "ENSEMBL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

head(converted_C9vsFUS)

reactome_genes_WTvsC9_entrez <- reactome_genes_WTvsC9[converted_WTvsC9$ENSEMBL]
names(reactome_genes_WTvsC9_entrez) <- converted_WTvsC9$ENTREZID
reactome_genes_WTvsC9_entrez <- sort(reactome_genes_WTvsC9_entrez, decreasing=TRUE)

reactome_genes_WTvsFUS_entrez <- reactome_genes_WTvsFUS[converted_WTvsFUS$ENSEMBL]
names(reactome_genes_WTvsFUS_entrez) <- converted_WTvsFUS$ENTREZID
reactome_genes_WTvsFUS_entrez <- sort(reactome_genes_WTvsFUS_entrez, decreasing=TRUE)

reactome_genes_C9vsFUS_entrez <- reactome_genes_C9vsFUS[converted_C9vsFUS$ENSEMBL]
names(reactome_genes_C9vsFUS_entrez) <- converted_C9vsFUS$ENTREZID
reactome_genes_C9vsFUS_entrez <- sort(reactome_genes_C9vsFUS_entrez, decreasing=TRUE)

reactome_WTvsC9 <- enrichPathway(gene = names(reactome_genes_WTvsC9_entrez), organism = 'human', pvalueCutoff = 0.05)
reactome_WTvsFUS <- enrichPathway(gene = names(reactome_genes_WTvsFUS_entrez), organism = 'human', pvalueCutoff = 0.05)
reactome_C9vsFUS <- enrichPathway(gene = names(reactome_genes_C9vsFUS_entrez), organism = 'human', pvalueCutoff = 0.05)

reactome_WTvsC9_sim <- pairwise_termsim(reactome_WTvsC9)
reactome_WTvsFUS_sim <- pairwise_termsim(reactome_WTvsFUS)
reactome_C9vsFUS_sim <- pairwise_termsim(reactome_C9vsFUS)

emapplot(reactome_WTvsC9_sim)
emapplot(reactome_WTvsFUS_sim)
emapplot(reactome_C9vsFUS_sim)


