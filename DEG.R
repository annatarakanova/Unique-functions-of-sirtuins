
######################### I. Differentially expression analysis ####

library(dplyr)
library(ggplot2)
library(DESeq2)
library(DOSE)
library(pheatmap)
library(tibble)
library(clusterProfiler)
library(rmarkdown)

# import raw feature count table 
data <- read.table("~/SIRT17_hMSCs/DESeq/salmon.merged.gene_counts.tsv", sep = "\t", header = TRUE)
head(data)

# filter of the low- and non-expressed genes (at least 10 counts per gene). Take only numeric columns 
data <- data[rowSums(data[, 3:ncol(data)]) > 10, ]

rownames(data) <- NULL

# create a meta table which maps our samples to the corresponding sample groups that we are investigating
# Extract sample names (excluding first two non-numeric columns)
sample_names <- colnames(data)[-(1:2)]
# Create metadata dataframe
meta <- data.frame(
  SampleID = ifelse(grepl("SIRT7", sample_names),
                    "SIRT7",
                    sub("_.*", "", sample_names)),  # For non-SIRT7 samples
  Genotype = ifelse(grepl("WT", sample_names), "WT", "KO"),
  Repeat = sub(".*_REP|.*_rep", "", sample_names),
  row.names = sample_names
)
# Convert Repeat to consistent uppercase format
meta$Repeat <- toupper(meta$Repeat)
# Reorder columns as requested
meta <- meta[, c("Genotype", "SampleID", "Repeat")]
# View the result
head(meta, n = 24)  # Shows all 24 samples

# Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(
  countData = round(data[, 3:ncol(data)]), 
  colData = meta,
  design = ~ Genotype  # Simple comparison of KO vs WT
)

dds


# Usually the step with estimateSizeFactors() is done automatically by DESeq2 during differential expression analysis.
# But let's do it by myself
# DESeq2 has a single estimateSizeFactors() function that will generate size factors for us)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
# retrieve the normalized counts matrix from dds
normalized_counts <- counts(dds, normalized=TRUE)


# Quality control
# the most simple method of log2 transformation
nts <- log2(assay(dds, normalized = TRUE)+1)
ggplot(as.data.frame(nts), aes(x = SIRT1_KO_REP1, y = SIRT1_KO_REP2)) + # plot two replicates 
  geom_point()+
  ggtitle('log2(x+1) transformation')+
  theme_dose(16)

# we see heteroscedasticity in the data -> perform rlog-transformation
rld <- rlog(dds, blind = TRUE)
ggplot(as.data.frame(assay(rld)), aes(x = SIRT1_KO_REP1, y = SIRT1_KO_REP2)) + 
  geom_point()+
  ggtitle('rlog transformation')+
  theme_dose(16)


install.packages("ggrepel") # package for labels 
library(ggrepel)
library(genefilter)

# PCA 
plotPCA.mystyle <- function (object, ntop = 500)
{
  font.size <- 18
  
  # --- 1. PCA Calculation (Unchanged) ---
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  
  # --- 2. Data Preparation for Plotting (Improved & Corrected) ---
  # Pulls metadata directly from the DESeq object for robustness
  plot_data <- data.frame(PC1 = pca$x[, 1], 
                          PC2 = pca$x[, 2], 
                          Genotype = colData(object)$Genotype, 
                          SampleID = colData(object)$SampleID)
  
  # --- 3. Create the ggplot (Updated) ---
  ggplot(data = plot_data, aes(x = PC1, y = PC2, color = Genotype)) +
    # Use geom_point to draw the points, colored by Genotype
    geom_point(size = 6) + 
    
    # Use geom_text_repel to add non-overlapping labels from the SampleID column
    geom_text_repel(aes(label = SampleID), 
                    box.padding = 0.5, # Increases space between label and point
                    max.overlaps = Inf) + # Ensures all labels are shown
    
    # Axis labels with variance explained
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    
    # Your custom theme elements
    # Note: I'm assuming theme_dose() is a custom function you have defined.
    # If not, you can replace it with a standard theme like theme_bw()
    # theme_dose(font.size = font.size) + 
    theme_bw(base_size = font.size) + # Using a standard theme for reproducibility
    theme(
      legend.key = element_rect(colour = NA, fill = NA), 
      legend.title = element_text(size = font.size - 2), # Title is useful
      legend.text = element_text(size = font.size - 2)
    )
}

# --- How to use the function ---
# Assuming 'rld' is your rlog-transformed object and its colData contains
# columns named "Genotype" and "SampleID".

# Plot PCA 
plotPCA.mystyle(rld)


# Correlation heatmap
# Extract the rlog matrix from the object
rld_mat <- assay(rld)
# Compute the pairwise correlation values for samples
rld_cor <- cor(rld_mat)
# Plot correlation values as a heatmap
pheatmap(rld_cor, annotation = meta[,c(1,2)])


# Differential testing 
meta$Factors <- paste(meta$Genotype, meta$SampleID, sep = '_') # join two columns with _ and create a new column Factors

# Create a new data type called factors. Factors are how R understands categorical data (i.e., groups).
# In statistical models like the one used by DESeq2, one of your groups must serve as the reference level or baseline. 
# All comparisons are made against this group.
meta$Factors <- factor(meta$Factors, levels = c("WT_WT", "KO_SIRT1", "KO_SIRT2", "KO_SIRT3", "KO_SIRT4", "KO_SIRT5", "KO_SIRT6", "KO_SIRT7"))

# Re-create the DESeqDataSet object
dds_factors <- DESeqDataSetFromMatrix(countData = round(data[, 3:ncol(data)]), 
                                      colData = meta, 
                                      design = ~ Factors)

# Check the model matrix in dds object
model.matrix(design(dds_factors), data = colData(dds_factors))

# DE analysis
dds_analysis <- DESeq(dds_factors)

# Plot dispersion estimates
plotDispEsts(dds_analysis)






######################## II. Find differentially expressed genes ####

# 1. Compare KO_SIRT1 with WT_WT
contrast_1 <- c("Factors", "KO_SIRT1", "WT_WT")
res_unshrunken_1 <- results(dds_analysis, contrast=contrast_1, alpha = 0.05)
res_1 <- lfcShrink(dds_analysis, 
                     contrast=contrast_1, 
                     res=res_unshrunken_1, 
                     type='normal')
knitr::kable(head(res_1))
plotMA(res_1, ylim=c(-3,5))
summary(res_1, alpha = 0.05)

# Extracting significant genes (|FDR p-value < 0.05| and |log2(Fold Change)| > 0.5)
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

# Convert the results table into a tibble
res_1_tb <- res_1 %>% 
  data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Select significant genes
sig_1 <- res_1_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
knitr::kable(head(sig_1))
nrow(sig_1) # 3162 significant genes


# Gene annotation
# In the initial table "data" we have a column with gene names. Let's connect gene numbers from
# sig_1 table with gene names from data table.
class(sig_1$gene) # output: character -> we need to convert them to numbers
class(res_1_tb$gene)
# Convert the 'gene' column to numeric
sig_1$gene <- as.numeric(as.character(sig_1$gene))
res_1_tb$gene <- as.numeric(as.character(res_1_tb$gene))
# Add a new column "gene_name" to sig_1
# by subsetting the "gene_name" column from "data" using the row numbers stored in sig_1$gene
sig_1$gene_name <- data$gene_name[sig_1$gene]
res_1_tb$gene_name <- data$gene_name[res_1_tb$gene]


# Volcano plot
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
p1 <- EnhancedVolcano(sig_1,
                      lab = sig_1$gene_name, 
                      x = 'log2FoldChange',
                      y = 'padj', 
                      title = 'KO vs WT, SIRT1',
                      subtitle = NULL,
                      pCutoff = 0.05, 
                      FCcutoff = 0.9)
cowplot::plot_grid(p1)





# 2. Compare KO_SIRT2 with WT_WT
contrast_2 <- c("Factors", "KO_SIRT2", "WT_WT")
res_unshrunken_2 <- results(dds_analysis, contrast=contrast_2, alpha = 0.05)
res_2 <- lfcShrink(dds_analysis, 
                   contrast=contrast_2, 
                   res=res_unshrunken_2, 
                   type='normal')
knitr::kable(head(res_2))
plotMA(res_2, ylim=c(-3,5))
summary(res_2, alpha = 0.05)

# Extracting significant genes (|FDR p-value < 0.05| and |log2(Fold Change)| > 0.5)
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

# Convert the results table into a tibble
res_2_tb <- res_2 %>% 
  data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Select significant genes
sig_2 <- res_2_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
knitr::kable(head(sig_2))
nrow(sig_2) # 3762 significant genes 


# Gene annotation
# In the initial table "data" we have a column with gene names. Let's connect gene numbers from
# sig_1 table with gene names from data table.
class(sig_2$gene) # output: character -> we need to convert them to numbers
class(res_2_tb$gene)
# Convert the 'gene' column to numeric
sig_2$gene <- as.numeric(as.character(sig_2$gene))
res_2_tb$gene <- as.numeric(as.character(res_2_tb$gene))
# Add a new column "gene_name" to sig_1
# by subsetting the "gene_name" column from "data" using the row numbers stored in sig_1$gene
sig_2$gene_name <- data$gene_name[sig_2$gene]
res_2_tb$gene_name <- data$gene_name[res_2_tb$gene]

# Volcano plot
# BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
p2 <- EnhancedVolcano(sig_2,
                      lab = sig_2$gene_name, 
                      x = 'log2FoldChange',
                      y = 'padj', 
                      title = 'KO vs WT, SIRT2',
                      subtitle = NULL,
                      pCutoff = 0.05, 
                      FCcutoff = 0.9)
cowplot::plot_grid(p2)





# 3. Compare KO_SIRT3 with WT_WT
contrast_3 <- c("Factors", "KO_SIRT3", "WT_WT")
res_unshrunken_3 <- results(dds_analysis, contrast=contrast_3, alpha = 0.05)
res_3 <- lfcShrink(dds_analysis, 
                   contrast=contrast_3, 
                   res=res_unshrunken_3, 
                   type='normal')
knitr::kable(head(res_3))
plotMA(res_3, ylim=c(-3,5))
summary(res_3, alpha = 0.05)

# Extracting significant genes (|FDR p-value < 0.05| and |log2(Fold Change)| > 0.5)
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

# Convert the results table into a tibble
res_3_tb <- res_3 %>% 
  data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Select significant genes
sig_3 <- res_3_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
knitr::kable(head(sig_3))
nrow(sig_3) # 2355 significant genes 


# Gene annotation
# In the initial table "data" we have a column with gene names. Let's connect gene numbers from
# sig_1 table with gene names from data table.
class(sig_3$gene) # output: character -> we need to convert them to numbers
class(res_3_tb$gene)
# Convert the 'gene' column to numeric
sig_3$gene <- as.numeric(as.character(sig_3$gene))
res_3_tb$gene <- as.numeric(as.character(res_3_tb$gene))
# Add a new column "gene_name" to sig_1
# by subsetting the "gene_name" column from "data" using the row numbers stored in sig_1$gene
sig_3$gene_name <- data$gene_name[sig_3$gene]
res_3_tb$gene_name <- data$gene_name[res_3_tb$gene]

# Volcano plot
# BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
p3 <- EnhancedVolcano(sig_3,
                      lab = sig_3$gene_name, 
                      x = 'log2FoldChange',
                      y = 'padj', 
                      title = 'KO vs WT, SIRT3',
                      subtitle = NULL,
                      pCutoff = 0.05, 
                      FCcutoff = 0.9)
cowplot::plot_grid(p3)





# 4. Compare KO_SIRT4 with WT_WT
contrast_4 <- c("Factors", "KO_SIRT4", "WT_WT")
res_unshrunken_4 <- results(dds_analysis, contrast=contrast_4, alpha = 0.05)
res_4 <- lfcShrink(dds_analysis, 
                   contrast=contrast_4, 
                   res=res_unshrunken_4, 
                   type='normal')
knitr::kable(head(res_4))
plotMA(res_4, ylim=c(-3,5))
summary(res_4, alpha = 0.05)

# Extracting significant genes (|FDR p-value < 0.05| and |log2(Fold Change)| > 0.5)
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

# Convert the results table into a tibble
res_4_tb <- res_4 %>% 
  data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Select significant genes
sig_4 <- res_4_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
knitr::kable(head(sig_4))
nrow(sig_4) # 2062 significant genes 


# Gene annotation
# In the initial table "data" we have a column with gene names. Let's connect gene numbers from
# sig_1 table with gene names from data table.
class(sig_4$gene) # output: character -> we need to convert them to numbers
class(res_4_tb$gene)
# Convert the 'gene' column to numeric
sig_4$gene <- as.numeric(as.character(sig_4$gene))
res_4_tb$gene <- as.numeric(as.character(res_4_tb$gene))
# Add a new column "gene_name" to sig_1
# by subsetting the "gene_name" column from "data" using the row numbers stored in sig_1$gene
sig_4$gene_name <- data$gene_name[sig_4$gene]
res_4_tb$gene_name <- data$gene_name[res_4_tb$gene]

# Volcano plot
# BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
p4 <- EnhancedVolcano(sig_4,
                      lab = sig_4$gene_name, 
                      x = 'log2FoldChange',
                      y = 'padj', 
                      title = 'KO vs WT, SIRT4',
                      subtitle = NULL,
                      pCutoff = 0.05, 
                      FCcutoff = 0.9)
cowplot::plot_grid(p4)





# 5. Compare KO_SIRT5 with WT_WT
contrast_5 <- c("Factors", "KO_SIRT5", "WT_WT")
res_unshrunken_5 <- results(dds_analysis, contrast=contrast_5, alpha = 0.05)
res_5 <- lfcShrink(dds_analysis, 
                   contrast=contrast_5, 
                   res=res_unshrunken_5, 
                   type='normal')
knitr::kable(head(res_5))
plotMA(res_5, ylim=c(-3,5))
summary(res_5, alpha = 0.05)

# Extracting significant genes (|FDR p-value < 0.05| and |log2(Fold Change)| > 0.5)
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

# Convert the results table into a tibble
res_5_tb <- res_5 %>% 
  data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Select significant genes
sig_5 <- res_5_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
knitr::kable(head(sig_5))
nrow(sig_5) # 2507 significant genes 


# Gene annotation
# In the initial table "data" we have a column with gene names. Let's connect gene numbers from
# sig_1 table with gene names from data table.
class(sig_5$gene) # output: character -> we need to convert them to numbers
class(res_5_tb$gene)
# Convert the 'gene' column to numeric
sig_5$gene <- as.numeric(as.character(sig_5$gene))
res_5_tb$gene <- as.numeric(as.character(res_5_tb$gene))
# Add a new column "gene_name" to sig_1
# by subsetting the "gene_name" column from "data" using the row numbers stored in sig_1$gene
sig_5$gene_name <- data$gene_name[sig_5$gene]
res_5_tb$gene_name <- data$gene_name[res_5_tb$gene]

# Volcano plot
# BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
p5 <- EnhancedVolcano(sig_5,
                      lab = sig_5$gene_name, 
                      x = 'log2FoldChange',
                      y = 'padj', 
                      title = 'KO vs WT, SIRT5',
                      subtitle = NULL,
                      pCutoff = 0.05, 
                      FCcutoff = 0.9)
cowplot::plot_grid(p5)





# 6. Compare KO_SIRT6 with WT_WT
contrast_6 <- c("Factors", "KO_SIRT6", "WT_WT")
res_unshrunken_6 <- results(dds_analysis, contrast=contrast_6, alpha = 0.05)
res_6 <- lfcShrink(dds_analysis, 
                   contrast=contrast_6, 
                   res=res_unshrunken_6, 
                   type='normal')
knitr::kable(head(res_6))
plotMA(res_6, ylim=c(-3,5))
summary(res_6, alpha = 0.05)

# Extracting significant genes (|FDR p-value < 0.05| and |log2(Fold Change)| > 0.5)
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

# Convert the results table into a tibble
res_6_tb <- res_6 %>% 
  data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Select significant genes
sig_6 <- res_6_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
knitr::kable(head(sig_6))
nrow(sig_6) # 3040 significant genes 


# Gene annotation
# In the initial table "data" we have a column with gene names. Let's connect gene numbers from
# sig_1 table with gene names from data table.
class(sig_6$gene) # output: character -> we need to convert them to numbers
class(res_6_tb$gene)
# Convert the 'gene' column to numeric
sig_6$gene <- as.numeric(as.character(sig_6$gene))
res_6_tb$gene <- as.numeric(as.character(res_6_tb$gene))
# Add a new column "gene_name" to sig_1
# by subsetting the "gene_name" column from "data" using the row numbers stored in sig_1$gene
sig_6$gene_name <- data$gene_name[sig_6$gene]
res_6_tb$gene_name <- data$gene_name[res_6_tb$gene]

# Volcano plot
# BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
p6 <- EnhancedVolcano(sig_6,
                      lab = sig_6$gene_name, 
                      x = 'log2FoldChange',
                      y = 'padj', 
                      title = 'KO vs WT, SIRT6',
                      subtitle = NULL,
                      pCutoff = 0.05, 
                      FCcutoff = 0.9)
cowplot::plot_grid(p6)





# 7. Compare KO_SIRT7 with WT_WT
contrast_7 <- c("Factors", "KO_SIRT7", "WT_WT")
res_unshrunken_7 <- results(dds_analysis, contrast=contrast_7, alpha = 0.05)
res_7 <- lfcShrink(dds_analysis, 
                   contrast=contrast_7, 
                   res=res_unshrunken_7, 
                   type='normal')
knitr::kable(head(res_7))
plotMA(res_7, ylim=c(-3,5))
summary(res_7, alpha = 0.05)

# Extracting significant genes (|FDR p-value < 0.05| and |log2(Fold Change)| > 0.5)
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

# Convert the results table into a tibble
res_7_tb <- res_7 %>% 
  data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Select significant genes
sig_7 <- res_7_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
knitr::kable(head(sig_7))
nrow(sig_7) # 2711 significant genes 


# Gene annotation
# In the initial table "data" we have a column with gene names. Let's connect gene numbers from
# sig_1 table with gene names from data table.
class(sig_7$gene) # output: character -> we need to convert them to numbers
class(res_7_tb$gene)
# Convert the 'gene' column to numeric
sig_7$gene <- as.numeric(as.character(sig_7$gene))
res_7_tb$gene <- as.numeric(as.character(res_7_tb$gene))
# Add a new column "gene_name" to sig_1
# by subsetting the "gene_name" column from "data" using the row numbers stored in sig_1$gene
sig_7$gene_name <- data$gene_name[sig_7$gene]
res_7_tb$gene_name <- data$gene_name[res_7_tb$gene]

# Volcano plot
# BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
p7 <- EnhancedVolcano(sig_7,
                      lab = sig_7$gene_name, 
                      x = 'log2FoldChange',
                      y = 'padj', 
                      title = 'KO vs WT, SIRT7',
                      subtitle = NULL,
                      pCutoff = 0.05, 
                      FCcutoff = 0.9)
cowplot::plot_grid(p7)







############################## IIIa. Find unique DEGs (gene names) ####

# Intersect sig dif genes from all the tables and find unique genes for each of the table sig_1, sig_2, ...
# Organizing data into a list
all_sig_tables <- list(
   SIRT1_KO = sig_1,
   SIRT2_KO = sig_2,
   SIRT3_KO = sig_3,
   SIRT4_KO = sig_4,
   SIRT5_KO = sig_5,
   SIRT6_KO = sig_6,
   SIRT7_KO = sig_7
 )

# It's cleaner to work with just the lists of gene names
all_gene_lists <- lapply(all_sig_tables, function(df) df$gene_name) # contain all the gene names from all the tables
print("Preview of the gene list for SIRT1_KO:")
print(head(all_gene_lists$SIRT1_KO))

# Loop through each list and find the unique genes
# Create an empty list to store the final results
unique_genes_list <- list()

# Loop through each element in our list of gene lists (from 1 to 7)
for (i in seq_along(all_gene_lists)) {
  
  # Get the name of the current condition (e.g., "SIRT1_KO")
  current_name <- names(all_gene_lists)[i]
  
  # STEP 1: pick one list
  current_genes <- all_gene_lists[[i]]
  
  # STEP 2: create a pool of genes of others lists 
  # The [-i] index is the R trick that means "select everything EXCEPT the i-th element"
  other_genes_pool <- unique(unlist(all_gene_lists[-i]))
  
  # STEP 3: compare gene list with pool of genes 
  # The setdiff() function finds all items in the first list that are NOT in the second list.
  unique_genes <- setdiff(current_genes, other_genes_pool)
  
  # Save the resulting vector of unique genes into our results list
  unique_genes_list[[current_name]] <- unique_genes
  
  # Print a message to show progress
  cat("Found", length(unique_genes), "unique genes for", current_name, "\n")
}

# View the number of unique genes found for each condition
summary_df <- data.frame(
  Condition = names(unique_genes_list),
  Unique_Gene_Count = sapply(unique_genes_list, length)
)
print(summary_df)





############################## IIIb. Find unique DEGs (gene tables) ####
# Create an empty list to store the final results WITH ALL COLUMNS
unique_genes_full_list <- list()

# Loop through each element in our list of gene lists (from 1 to 7)
for (i in seq_along(all_sig_tables)) {
  
  # Get the name of the current condition (e.g., "SIRT1_KO")
  current_name <- names(all_sig_tables)[i]
  
  # STEP 1: Get the current table with all columns
  current_table <- all_sig_tables[[i]]
  
  # STEP 2: Create a pool of genes from other lists 
  other_genes_pool <- unique(unlist(all_gene_lists[-i]))
  
  # STEP 3: Filter the current table to keep only unique genes
  # Using the gene_name column to identify unique genes
  unique_genes_full <- current_table %>% 
    filter(!gene_name %in% other_genes_pool)
  
  # Save the resulting data frame with all columns into our results list
  unique_genes_full_list[[current_name]] <- unique_genes_full
  
  # Print a message to show progress
  cat("Found", nrow(unique_genes_full), "unique genes for", current_name, "\n")
}

# View the number of unique genes found for each condition
summary_df <- data.frame(
  Condition = names(unique_genes_full_list),
  Unique_Gene_Count = sapply(unique_genes_full_list, nrow)
)
print(summary_df)

# You can now access the full unique gene tables for each SIRT:
# For example, to see unique genes for SIRT3 with all columns:
print("Unique genes for SIRT3_KO with all columns:")
print(head(unique_genes_full_list$SIRT3_KO))





########## IV. GO, KEGG, Reactome, Msigdbr enrichment analysis on unique DEGs ####

# Functional analysis
library(dplyr)
# Gene ontology (GO) analysis 
library(clusterProfiler)
library(org.Hs.eg.db) #  Homo sapiens (Hs) database

# Use all the expressed genes as a background 
background <- unique(c(
  res_1_tb$gene_name, res_2_tb$gene_name, res_3_tb$gene_name,
  res_4_tb$gene_name, res_5_tb$gene_name, res_6_tb$gene_name,
  res_7_tb$gene_name
))


# 1. SIRT1_KO

# Get the list of significant genes for SIRT1_KO
sirt1_genes <- sig_1$gene_name

# Unique genes for SIRT1_KO
unique_to_sirt1 <- unique_genes_list$SIRT1_KO
# Define the name for your output text file
output_filename <- "unique_sirt1_genes.txt"
output_directory <- "~/SIRT17_hMSCs/DESeq"
full_path_to_file <- file.path(output_directory, output_filename)
# Save the file in current working directory
write(
  unique_to_sirt1,              # The vector of gene names to write
  file = full_path_to_file,       # The name of the file to create
  sep = "\n"                    # The separator to use between elements
)

# Let's check how many we found
cat("Total significant genes in SIRT1_KO:", length(sirt1_genes), "\n")
cat("Genes unique to SIRT1_KO:", length(unique_to_sirt1), "\n")

# Perform the GO analysis for the SIRT1_KO significant genes
ego_1 <- enrichGO(
  gene          = unique_to_sirt1,
  universe      = background,
  keyType       = "SYMBOL",          # We are using gene symbols
  OrgDb         = org.Hs.eg.db,      # Using the Human organism database
  ont           = "BP",              # Analyzing Biological Processes
  pAdjustMethod = "BH",              # Benjamini-Hochberg for p-value correction
  pvalueCutoff  = 1,
  qvalueCutoff = 0.5,
  readable      = TRUE               # This is a useful extra argument!
)

# d.ego_1 <- data.frame(ego_1@result)

# View the enriched categories
p1 <- barplot(ego_1, 
              # showCategory = 10, 
              title = "GO analysis for unique SIRT1-KO genes")

# Display the plot
cowplot::plot_grid(p1) # X-axes - number of genes that are associated with each biological process


# Pathway analysis (KEGG)
# --- Step 1: Convert the 'gene' list from SYMBOL to ENTREZID ---

# Use bitr to translate the symbols in the unique_to_sirt1 vector
# bitr returns a data frame with columns for the original and translated IDs
sirt1_ids_df <- bitr(unique_to_sirt1,          # Your vector of gene symbols
                     fromType = "SYMBOL",      # The type of ID you are providing
                     toType = "ENTREZID",      # The type of ID you want to get
                     OrgDb = org.Hs.eg.db)

# Now, extract the ENTREZID column from this new data frame
sirt1_entrez_ids <- sirt1_ids_df$ENTREZID


# --- Step 2: Convert the 'universe' list from SYMBOL to ENTREZID ---

# Repeat the process for your universe gene list
background_ids_df <- bitr(background,
                        fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)

# Extract the ENTREZID column
background_entrez_ids <- background_ids_df$ENTREZID


# --- Step 3: Run enrichKEGG with the CORRECTED Entrez ID vectors ---
# Now your inputs are simple vectors of Entrez IDs, and na.omit() will work correctly.

ekegg_1 <- enrichKEGG(gene = na.omit(sirt1_entrez_ids),
                      organism = 'hsa',
                      pvalueCutoff = 1,
                      qvalueCutoff = 0.5,
                      universe = na.omit(background_entrez_ids),
                      minGSSize = 2 ) # each pathway should be enriched by min 2 genes
                

# ekegg_readable <- setReadable(
  # ekegg_1,               # Your enrichKEGG result object
  # OrgDb = org.Hs.eg.db,        # The organism database for mapping
  # keyType = "ENTREZID"         # The type of ID in your object
# )

# Visualize the results 
library(enrichplot)
ups_1 <- enrichplot::upsetplot(ekegg_1)
ups_1



# Reactome analysis
BiocManager::install("ReactomePA")
library(ReactomePA)

# The input gene ID should be Entrez gene ID
pathway_results <- enrichPathway(
  gene          = sirt1_entrez_ids,
  organism      = "human",
  universe      = background_entrez_ids,
  pvalueCutoff  = 1,
  qvalueCutoff = 0.5,
  pAdjustMethod = "BH",
  readable      = TRUE
)

# head(pathway_results)

# View the enriched categories
p1 <- barplot(pathway_results, 
              # showCategory = 15, 
              title = "Reactome analysis for unique SIRT1-KO genes")

# Display the plot
cowplot::plot_grid(p1)



# msigdbr - signature enrichment analysis (можно посмотреть не только онтологии, но и разные сигнатуры)
install.packages("msigdbr")
library(msigdbr)

# Шаг 1: Анализ с коллекцией Hallmark (H) - Общий обзор
# 1. Загружаем коллекцию Hallmark
hallmark_gsets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

# 2. Запускаем анализ обогащения
enrich_hallmark <- enricher(
  gene = sirt1_entrez_ids,
  universe = background_entrez_ids,
  pvalueCutoff = 1,
  qvalueCutoff = 0.5,
  TERM2GENE = hallmark_gsets
)

# 3. Визуализируем результаты (если они есть)
if (!is.null(enrich_hallmark) && nrow(enrich_hallmark) > 0) {
  p_hallmark <- dotplot(enrich_hallmark, showCategory = 15) +
    ggtitle("Hallmark Pathway Enrichment, SIRT1")
  print(p_hallmark)
} else {
  cat("Не найдено значимых путей в коллекции Hallmark.\n")
}



# Шаг 2: Углубленный анализ путей с C2:CP:REACTOME
# 1. Загружаем коллекцию путей Reactome
reactome_gsets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, entrez_gene)

# 2. Запускаем анализ
enrich_reactome <- enricher(
  gene = sirt1_entrez_ids,
  universe = background_entrez_ids,
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  TERM2GENE = reactome_gsets
)

# 3. Визуализируем
if (!is.null(enrich_reactome) && nrow(enrich_reactome) > 0) {
  p_reactome <- dotplot(enrich_reactome, showCategory = 15) +
    ggtitle("Reactome Pathway Enrichment, SIRT1")
  print(p_reactome)
} else {
  cat("Не найдено значимых путей в коллекции Reactome.\n")
}

p_reactome <- dotplot(
  enrich_reactome, 
  showCategory = 15, 
  title = "Reactome Pathway Enrichment, SIRT1"
)

# --- 2. Добавьте слой theme() для изменения размера шрифта ---
# Мы изменяем размер шрифта для текста на оси Y (axis.text.y)
# и на оси X (axis.text.x) для лучшего вида.

p_reactome_final <- p_reactome + 
  theme(
    axis.text.y = element_text(size = 4),  # Устанавливаем размер шрифта для названий путей
    axis.text.x = element_text(size = 10)  # Можно также немного уменьшить шрифт на оси X
  )

# --- 3. Отобразите итоговый график ---
print(p_reactome_final)



# Шаг 3: Поиск регуляторов с C3:TFT (Мишени транскрипционных факторов)
# 1. Загружаем коллекцию мишеней ТФ
tft_gsets <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD") %>%
  dplyr::select(gs_name, entrez_gene)

# 2. Запускаем анализ
enrich_tft <- enricher(
  gene = sirt1_entrez_ids,
  universe = background_entrez_ids,
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  TERM2GENE = tft_gsets
)

# 3. Визуализируем
if (!is.null(enrich_tft) && nrow(enrich_tft) > 0) {
  p_tft <- dotplot(enrich_tft, showCategory = 15) +
    ggtitle("Transcription Factor Target Enrichment, SIRT1")
  print(p_tft)
} else {
  cat("Не найдено значимого обогащения по мишеням транскрипционных факторов.\n")
}



# Шаг 4: Стандартный GO-анализ с C5:GO:BP (Биологические процессы)
# 1. Загружаем коллекцию GO:BP
go_bp_gsets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  dplyr::select(gs_name, entrez_gene)

# 2. Запускаем анализ
enrich_go_bp <- enricher(
  gene = sirt1_entrez_ids,
  universe = background_entrez_ids,
  pvalueCutoff = 0.05,
  qvalueCutoff = 1,
  TERM2GENE = go_bp_gsets
)

# 3. Визуализируем
if (!is.null(enrich_go_bp) && nrow(enrich_go_bp) > 0) {
  # Можно применить simplify() для удаления избыточных GO-терминов
  enrich_go_bp_simple <- simplify(enrich_go_bp, cutoff = 0.7, by = "p.adjust")
  
  p_go_bp <- dotplot(enrich_go_bp_simple, showCategory = 15) +
    ggtitle("Gene Ontology (BP) Enrichment")
  print(p_go_bp)
} else {
  cat("Не найдено значимых путей в коллекции Gene Ontology (BP).\n")
}






# 2. SIRT2_KO

# Get the list of significant genes for SIRT2_KO
sirt2_genes <- sig_2$gene_name

# Unique genes for SIRT2_KO
unique_to_sirt2 <- unique_genes_list$SIRT2_KO
# Define the name for your output text file
output_filename <- "unique_sirt2_genes.txt"
output_directory <- "~/SIRT17_hMSCs/DESeq"
full_path_to_file <- file.path(output_directory, output_filename)
# Save the file in current working directory
write(
  unique_to_sirt2,              # The vector of gene names to write
  file = full_path_to_file,       # The name of the file to create
  sep = "\n"                    # The separator to use between elements
)

# Let's check how many we found
cat("Total significant genes in SIRT2_KO:", length(sirt2_genes), "\n")
cat("Genes unique to SIRT2_KO:", length(unique_to_sirt2), "\n")

# Perform the GO analysis for the SIRT2_KO significant genes
ego_2 <- enrichGO(
  gene          = unique_to_sirt2,
  universe      = background,
  keyType       = "SYMBOL",          # We are using gene symbols
  OrgDb         = org.Hs.eg.db,      # Using the Human organism database
  ont           = "BP",              # Analyzing Biological Processes
  pAdjustMethod = "BH",              # Benjamini-Hochberg for p-value correction
  pvalueCutoff  = 1,
  qvalueCutoff = 1,
  readable      = TRUE               # This is a useful extra argument!
)

# View the enriched categories
p2 <- barplot(ego_2, title = "GO analysis for unique SIRT2-KO genes")

# Display the plot
cowplot::plot_grid(p2)


# Pathway analysis
# --- Step 1: Convert the 'gene' list from SYMBOL to ENTREZID ---

# Use bitr to translate the symbols in the unique_to_sirt2 vector
# bitr returns a data frame with columns for the original and translated IDs
sirt2_ids_df <- bitr(unique_to_sirt2,          # Your vector of gene symbols
                     fromType = "SYMBOL",      # The type of ID you are providing
                     toType = "ENTREZID",      # The type of ID you want to get
                     OrgDb = org.Hs.eg.db)

# Now, extract the ENTREZID column from this new data frame
sirt2_entrez_ids <- sirt2_ids_df$ENTREZID


# --- Step 2: Convert the 'universe' list from SYMBOL to ENTREZID ---

# Repeat the process for your universe gene list
background_ids_df <- bitr(background,
                        fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)

# Extract the ENTREZID column
background_entrez_ids <- background_ids_df$ENTREZID


# --- Step 3: Run enrichKEGG with the CORRECTED Entrez ID vectors ---
# Now your inputs are simple vectors of Entrez IDs, and na.omit() will work correctly.

ekegg_2 <- enrichKEGG(gene = na.omit(sirt2_entrez_ids),
                      organism = 'hsa',
                      pvalueCutoff = 1,
                      qvalueCutoff = 0.5,
                      universe = na.omit(background_entrez_ids))


# View the results (if any are found)
head(ekegg_2)

library(enrichplot)
ups_2 <- enrichplot::upsetplot(ekegg_2)
ups_2


# Reactome analysis
BiocManager::install("ReactomePA")
library(ReactomePA)

# The input gene ID should be Entrez gene ID
pathway_results <- enrichPathway(
  gene          = sirt2_entrez_ids,
  organism      = "human",
  universe      = background_entrez_ids,
  pvalueCutoff  = 1,
  qvalueCutoff  = 0.5,
  pAdjustMethod = "BH",
  readable      = TRUE
)

head(pathway_results)

# View the enriched categories
p2 <- barplot(pathway_results, 
              # showCategory = 15, 
              title = "Reactome analysis for unique SIRT2-KO genes")

# Display the plot
cowplot::plot_grid(p2)


# msigdbr - signature enrichment analysis (можно посмотреть не только онтологии, но и разные сигнатуры)
# Шаг 1: Анализ с коллекцией Hallmark (H) - Общий обзор
# 1. Загружаем коллекцию Hallmark
hallmark_gsets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

# 2. Запускаем анализ обогащения
enrich_hallmark <- enricher(
  gene = sirt2_entrez_ids,
  universe = background_entrez_ids,
  pvalueCutoff = 1,
  qvalueCutoff = 0.5,
  TERM2GENE = hallmark_gsets
)

# 3. Визуализируем результаты (если они есть)
if (!is.null(enrich_hallmark) && nrow(enrich_hallmark) > 0) {
  p_hallmark <- dotplot(enrich_hallmark, showCategory = 15) +
    ggtitle("Hallmark Pathway Enrichment, SIRT2")
  print(p_hallmark)
} else {
  cat("Не найдено значимых путей в коллекции Hallmark.\n")
}




# 3. SIRT3_KO

# Get the list of significant genes for SIRT3_KO
sirt3_genes <- sig_3$gene_name

# Unique genes for SIRT3_KO
unique_to_sirt3 <- unique_genes_list$SIRT3_KO
# Define the name for your output text file
output_filename <- "unique_sirt3_genes.txt"
output_directory <- "~/SIRT17_hMSCs/DESeq"
full_path_to_file <- file.path(output_directory, output_filename)
# Save the file in current working directory
write(
  unique_to_sirt3,              # The vector of gene names to write
  file = full_path_to_file,       # The name of the file to create
  sep = "\n"                    # The separator to use between elements
)

# Let's check how many we found
cat("Total significant genes in SIRT3_KO:", length(sirt3_genes), "\n")
cat("Genes unique to SIRT3_KO:", length(unique_to_sirt3), "\n")

# Perform the GO analysis for the SIRT3_KO significant genes
ego_3 <- enrichGO(
  gene          = unique_to_sirt3,
  universe      = background,
  keyType       = "SYMBOL",          # We are using gene symbols
  OrgDb         = org.Hs.eg.db,      # Using the Human organism database
  ont           = "BP",              # Analyzing Biological Processes
  pAdjustMethod = "BH",              # Benjamini-Hochberg for p-value correction
  pvalueCutoff  = 1,
  qvalueCutoff  = 0.5,
  readable      = TRUE               # This is a useful extra argument!
)

# View the enriched categories
p3 <- barplot(ego_3,
              showCategory = 6,
              title = "GO analysis for unique SIRT3-KO genes")

# Display the plot
cowplot::plot_grid(p3)


# Pathway analysis
# --- Step 1: Convert the 'gene' list from SYMBOL to ENTREZID ---

# Use bitr to translate the symbols in the unique_to_sirt2 vector
# bitr returns a data frame with columns for the original and translated IDs
sirt3_ids_df <- bitr(unique_to_sirt3,          # Your vector of gene symbols
                     fromType = "SYMBOL",      # The type of ID you are providing
                     toType = "ENTREZID",      # The type of ID you want to get
                     OrgDb = org.Hs.eg.db)

# Now, extract the ENTREZID column from this new data frame
sirt3_entrez_ids <- sirt3_ids_df$ENTREZID


# --- Step 2: Convert the 'universe' list from SYMBOL to ENTREZID ---

# Repeat the process for your universe gene list
background_ids_df <- bitr(background,
                          fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)

# Extract the ENTREZID column
background_entrez_ids <- background_ids_df$ENTREZID


# --- Step 3: Run enrichKEGG with the CORRECTED Entrez ID vectors ---
# Now your inputs are simple vectors of Entrez IDs, and na.omit() will work correctly.

ekegg_3 <- enrichKEGG(gene = na.omit(sirt3_entrez_ids),
                      organism = 'hsa',
                      pvalueCutoff = 1,
                      qvalueCutoff = 0.5,
                      universe = na.omit(background_entrez_ids))


# View the results (if any are found)
head(ekegg_3)

library(enrichplot)
ups_3 <- enrichplot::upsetplot(ekegg_3)
ups_3


# Reactome analysis
BiocManager::install("ReactomePA")
library(ReactomePA)

# The input gene ID should be Entrez gene ID
pathway_results <- enrichPathway(
  gene          = sirt3_entrez_ids,
  organism      = "human",
  universe      = background_entrez_ids,
  pvalueCutoff  = 1,
  qvalueCutoff  = 0.5,
  pAdjustMethod = "BH",
  readable      = TRUE
)

head(pathway_results)

# View the enriched categories
p3 <- barplot(pathway_results, 
              # showCategory = 15, 
              title = "Reactome analysis for unique SIRT3-KO genes")

# Display the plot
cowplot::plot_grid(p3)



# msigdbr - signature enrichment analysis (можно посмотреть не только онтологии, но и разные сигнатуры)
# Шаг 1: Анализ с коллекцией Hallmark (H) - Общий обзор
# 1. Загружаем коллекцию Hallmark
hallmark_gsets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

# 2. Запускаем анализ обогащения
enrich_hallmark <- enricher(
  gene = sirt3_entrez_ids,
  universe = background_entrez_ids,
  pvalueCutoff = 1,
  qvalueCutoff = 0.5,
  TERM2GENE = hallmark_gsets
)

# 3. Визуализируем результаты (если они есть)
if (!is.null(enrich_hallmark) && nrow(enrich_hallmark) > 0) {
  p_hallmark <- dotplot(enrich_hallmark, showCategory = 15) +
    ggtitle("Hallmark Pathway Enrichment, SIRT3")
  print(p_hallmark)
} else {
  cat("Не найдено значимых путей в коллекции Hallmark.\n")
}





# 4. SIRT4_KO

# Get the list of significant genes for SIRT4_KO
sirt4_genes <- sig_4$gene_name

# Unique genes for SIRT4_KO
unique_to_sirt4 <- unique_genes_list$SIRT4_KO
# Define the name for your output text file
output_filename <- "unique_sirt4_genes.txt"
output_directory <- "~/SIRT17_hMSCs/DESeq"
full_path_to_file <- file.path(output_directory, output_filename)
# Save the file in current working directory
write(
  unique_to_sirt4,              # The vector of gene names to write
  file = full_path_to_file,       # The name of the file to create
  sep = "\n"                    # The separator to use between elements
)

# Let's check how many we found
cat("Total significant genes in SIRT4_KO:", length(sirt4_genes), "\n")
cat("Genes unique to SIRT4_KO:", length(unique_to_sirt4), "\n")

# Perform the GO analysis for the SIRT4_KO significant genes
ego_4 <- enrichGO(
  gene          = unique_to_sirt4,
  universe      = background,
  keyType       = "SYMBOL",          # We are using gene symbols
  OrgDb         = org.Hs.eg.db,      # Using the Human organism database
  ont           = "BP",              # Analyzing Biological Processes
  pAdjustMethod = "BH",              # Benjamini-Hochberg for p-value correction
  pvalueCutoff  = 1,
  qvalueCutoff  = 0.5,
  readable      = TRUE               # This is a useful extra argument!
)

# View the enriched categories
p4 <- barplot(ego_4, title = "GO analysis for unique SIRT4-KO genes")

# Display the plot
cowplot::plot_grid(p4)



# Convert the 'gene' list from SYMBOL to ENTREZID

# Use bitr to translate the symbols in the unique_to_sirt4 vector
# bitr returns a data frame with columns for the original and translated IDs
sirt4_ids_df <- bitr(unique_to_sirt4,          # Your vector of gene symbols
                     fromType = "SYMBOL",      # The type of ID you are providing
                     toType = "ENTREZID",      # The type of ID you want to get
                     OrgDb = org.Hs.eg.db)

# Now, extract the ENTREZID column from this new data frame
sirt4_entrez_ids <- sirt4_ids_df$ENTREZID

# --- Step 3: Run enrichKEGG with the CORRECTED Entrez ID vectors ---
# Now your inputs are simple vectors of Entrez IDs, and na.omit() will work correctly.

ekegg_4 <- enrichKEGG(gene = na.omit(sirt4_entrez_ids),
                      organism = 'hsa',
                      pvalueCutoff = 1,
                      qvalueCutoff = 0.5,
                      universe = na.omit(background_entrez_ids))


# View the results (if any are found)
head(ekegg_4)

library(enrichplot)
ups_4 <- enrichplot::upsetplot(ekegg_4)
ups_4


# Reactome analysis
BiocManager::install("ReactomePA")
library(ReactomePA)

# The input gene ID should be Entrez gene ID
pathway_results <- enrichPathway(
  gene          = sirt4_entrez_ids,
  organism      = "human",
  universe      = background_entrez_ids,
  pvalueCutoff  = 1,
  qvalueCutoff  = 0.5,
  pAdjustMethod = "BH",
  readable      = TRUE
)

head(pathway_results)

# View the enriched categories
p4 <- barplot(pathway_results, 
              # showCategory = 15, 
              title = "Reactome analysis for unique SIRT4-KO genes")

# Display the plot
cowplot::plot_grid(p4)



# msigdbr - signature enrichment analysis (можно посмотреть не только онтологии, но и разные сигнатуры)
# Шаг 1: Анализ с коллекцией Hallmark (H) - Общий обзор
# 1. Загружаем коллекцию Hallmark
hallmark_gsets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

# 2. Запускаем анализ обогащения
enrich_hallmark <- enricher(
  gene = sirt4_entrez_ids,
  universe = background_entrez_ids,
  pvalueCutoff = 1,
  qvalueCutoff = 0.5,
  TERM2GENE = hallmark_gsets
)

# 3. Визуализируем результаты (если они есть)
if (!is.null(enrich_hallmark) && nrow(enrich_hallmark) > 0) {
  p_hallmark <- dotplot(enrich_hallmark, showCategory = 15) +
    ggtitle("Hallmark Pathway Enrichment, SIRT4")
  print(p_hallmark)
} else {
  cat("Не найдено значимых путей в коллекции Hallmark.\n")
}





# 5. SIRT5_KO

# Get the list of significant genes for SIRT5_KO
sirt5_genes <- sig_5$gene_name

# Unique genes for SIRT5_KO
unique_to_sirt5 <- unique_genes_list$SIRT5_KO
# Define the name for your output text file
output_filename <- "unique_sirt5_genes.txt"
output_directory <- "~/SIRT17_hMSCs/DESeq"
full_path_to_file <- file.path(output_directory, output_filename)
# Save the file in current working directory
write(
  unique_to_sirt5,              # The vector of gene names to write
  file = full_path_to_file,       # The name of the file to create
  sep = "\n"                    # The separator to use between elements
)

# Let's check how many we found
cat("Total significant genes in SIRT5_KO:", length(sirt5_genes), "\n")
cat("Genes unique to SIRT5_KO:", length(unique_to_sirt5), "\n")

# Perform the GO analysis for the SIRT5_KO significant genes
ego_5 <- enrichGO(
  gene          = unique_to_sirt5,
  universe      = background,
  keyType       = "SYMBOL",          # We are using gene symbols
  OrgDb         = org.Hs.eg.db,      # Using the Human organism database
  ont           = "BP",              # Analyzing Biological Processes
  pAdjustMethod = "BH",              # Benjamini-Hochberg for p-value correction
  pvalueCutoff  = 1,
  qvalueCutoff  = 0.5,
  readable      = TRUE               # This is a useful extra argument!
)

# View the enriched categories
p5 <- barplot(ego_5, title = "GO analysis for unique SIRT5-KO genes")

# Display the plot
cowplot::plot_grid(p5)


# Convert the 'gene' list from SYMBOL to ENTREZID

# Use bitr to translate the symbols in the unique_to_sirt5 vector
# bitr returns a data frame with columns for the original and translated IDs
sirt5_ids_df <- bitr(unique_to_sirt5,          # Your vector of gene symbols
                     fromType = "SYMBOL",      # The type of ID you are providing
                     toType = "ENTREZID",      # The type of ID you want to get
                     OrgDb = org.Hs.eg.db)

# Now, extract the ENTREZID column from this new data frame
sirt5_entrez_ids <- sirt5_ids_df$ENTREZID

# --- Step 3: Run enrichKEGG with the CORRECTED Entrez ID vectors ---
# Now your inputs are simple vectors of Entrez IDs, and na.omit() will work correctly.

ekegg_5 <- enrichKEGG(gene = na.omit(sirt5_entrez_ids),
                      organism = 'hsa',
                      pvalueCutoff = 1,
                      qvalueCutoff = 0.5,
                      universe = na.omit(background_entrez_ids))


# View the results (if any are found)
head(ekegg_5)

library(enrichplot)
ups_5 <- enrichplot::upsetplot(ekegg_5)
ups_5


# Reactome analysis
BiocManager::install("ReactomePA")
library(ReactomePA)

# The input gene ID should be Entrez gene ID
pathway_results <- enrichPathway(
  gene          = sirt5_entrez_ids,
  organism      = "human",
  universe      = background_entrez_ids,
  pvalueCutoff  = 1,
  qvalueCutoff  = 0.5,
  pAdjustMethod = "BH",
  readable      = TRUE
)

head(pathway_results)

# View the enriched categories
p5 <- barplot(pathway_results, 
              showCategory = 7, 
              title = "Reactome analysis for unique SIRT5-KO genes")

# Display the plot
cowplot::plot_grid(p5)





# msigdbr - signature enrichment analysis (можно посмотреть не только онтологии, но и разные сигнатуры)
# Шаг 1: Анализ с коллекцией Hallmark (H) - Общий обзор
# 1. Загружаем коллекцию Hallmark
hallmark_gsets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

# 2. Запускаем анализ обогащения
enrich_hallmark <- enricher(
  gene = sirt5_entrez_ids,
  universe = background_entrez_ids,
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  TERM2GENE = hallmark_gsets
)

# 3. Визуализируем результаты (если они есть)
if (!is.null(enrich_hallmark) && nrow(enrich_hallmark) > 0) {
  p_hallmark <- dotplot(enrich_hallmark, showCategory = 15) +
    ggtitle("Hallmark Pathway Enrichment, SIRT5")
  print(p_hallmark)
} else {
  cat("Не найдено значимых путей в коллекции Hallmark.\n")
}





# 6. SIRT6_KO

# Get the list of significant genes for SIRT6_KO
sirt6_genes <- sig_6$gene_name

# Unique genes for SIRT6_KO
unique_to_sirt6 <- unique_genes_list$SIRT6_KO
# Define the name for your output text file
output_filename <- "unique_sirt6_genes.txt"
output_directory <- "~/SIRT17_hMSCs/DESeq"
full_path_to_file <- file.path(output_directory, output_filename)
# Save the file in current working directory
write(
  unique_to_sirt6,              # The vector of gene names to write
  file = full_path_to_file,       # The name of the file to create
  sep = "\n"                    # The separator to use between elements
)

# Let's check how many we found
cat("Total significant genes in SIRT6_KO:", length(sirt6_genes), "\n")
cat("Genes unique to SIRT6_KO:", length(unique_to_sirt6), "\n")

# Perform the GO analysis for the SIRT6_KO significant genes
ego_6 <- enrichGO(
  gene          = unique_to_sirt6,
  universe      = background,
  keyType       = "SYMBOL",          # We are using gene symbols
  OrgDb         = org.Hs.eg.db,      # Using the Human organism database
  ont           = "BP",              # Analyzing Biological Processes
  pAdjustMethod = "BH",              # Benjamini-Hochberg for p-value correction
  pvalueCutoff  = 1,
  qvalueCutoff  = 0.5,
  readable      = TRUE               # This is a useful extra argument!
)

# View the enriched categories
p6 <- barplot(ego_6, title = "GO analysis for unique SIRT6-KO genes")

# Display the plot
cowplot::plot_grid(p6)


# Convert the 'gene' list from SYMBOL to ENTREZID

# Use bitr to translate the symbols in the unique_to_sirt6 vector
# bitr returns a data frame with columns for the original and translated IDs
sirt6_ids_df <- bitr(unique_to_sirt6,          # Your vector of gene symbols
                     fromType = "SYMBOL",      # The type of ID you are providing
                     toType = "ENTREZID",      # The type of ID you want to get
                     OrgDb = org.Hs.eg.db)

# Now, extract the ENTREZID column from this new data frame
sirt6_entrez_ids <- sirt6_ids_df$ENTREZID

# --- Step 3: Run enrichKEGG with the CORRECTED Entrez ID vectors ---
# Now your inputs are simple vectors of Entrez IDs, and na.omit() will work correctly.

ekegg_6 <- enrichKEGG(gene = na.omit(sirt6_entrez_ids),
                      organism = 'hsa',
                      pvalueCutoff = 1,
                      qvalueCutoff = 0.5,
                      universe = na.omit(background_entrez_ids))


# View the results (if any are found)
head(ekegg_6)

library(enrichplot)
ups_6 <- enrichplot::upsetplot(ekegg_6)
ups_6


# Reactome analysis
BiocManager::install("ReactomePA")
library(ReactomePA)

# The input gene ID should be Entrez gene ID
pathway_results <- enrichPathway(
  gene          = sirt6_entrez_ids,
  organism      = "human",
  universe      = background_entrez_ids,
  pvalueCutoff  = 1,
  qvalueCutoff  = 0.5,
  pAdjustMethod = "BH",
  readable      = TRUE
)

head(pathway_results)

# View the enriched categories
p6 <- barplot(pathway_results, 
              # showCategory = 7, 
              title = "Reactome analysis for unique SIRT6-KO genes")

# Display the plot
cowplot::plot_grid(p6)






# msigdbr - signature enrichment analysis (можно посмотреть не только онтологии, но и разные сигнатуры)
# Шаг 1: Анализ с коллекцией Hallmark (H) - Общий обзор
# 1. Загружаем коллекцию Hallmark
hallmark_gsets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

# 2. Запускаем анализ обогащения
enrich_hallmark <- enricher(
  gene = sirt6_entrez_ids,
  universe = background_entrez_ids,
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  TERM2GENE = hallmark_gsets
)

# 3. Визуализируем результаты (если они есть)
if (!is.null(enrich_hallmark) && nrow(enrich_hallmark) > 0) {
  p_hallmark <- dotplot(enrich_hallmark, showCategory = 15) +
    ggtitle("Hallmark Pathway Enrichment, SIRT6")
  print(p_hallmark)
} else {
  cat("Не найдено значимых путей в коллекции Hallmark.\n")
}






# 7. SIRT7_KO

# Get the list of significant genes for SIRT7_KO
sirt7_genes <- sig_7$gene_name

# Unique genes for SIRT7_KO
unique_to_sirt7 <- unique_genes_list$SIRT7_KO
# Define the name for your output text file
output_filename <- "unique_sirt7_genes.txt"
output_directory <- "~/SIRT17_hMSCs/DESeq"
full_path_to_file <- file.path(output_directory, output_filename)
# Save the file in current working directory
write(
  unique_to_sirt7,              # The vector of gene names to write
  file = full_path_to_file,       # The name of the file to create
  sep = "\n"                    # The separator to use between elements
)

# Let's check how many we found
cat("Total significant genes in SIRT7_KO:", length(sirt7_genes), "\n")
cat("Genes unique to SIRT7_KO:", length(unique_to_sirt7), "\n")

# Perform the GO analysis for the SIRT7_KO significant genes
ego_7 <- enrichGO(
  gene          = unique_to_sirt7,
  universe      = background,
  keyType       = "SYMBOL",          # We are using gene symbols
  OrgDb         = org.Hs.eg.db,      # Using the Human organism database
  ont           = "BP",              # Analyzing Biological Processes
  pAdjustMethod = "BH",              # Benjamini-Hochberg for p-value correction
  pvalueCutoff  = 1,
  qvalueCutoff  = 0.5,
  readable      = TRUE               # This is a useful extra argument!
)

# View the enriched categories
p7 <- barplot(ego_7, title = "GO analysis for unique SIRT7-KO genes")

# Display the plot
cowplot::plot_grid(p7)



#Convert the 'gene' list from SYMBOL to ENTREZID

# Use bitr to translate the symbols in the unique_to_sirt7 vector
# bitr returns a data frame with columns for the original and translated IDs
sirt7_ids_df <- bitr(unique_to_sirt7,          # Your vector of gene symbols
                     fromType = "SYMBOL",      # The type of ID you are providing
                     toType = "ENTREZID",      # The type of ID you want to get
                     OrgDb = org.Hs.eg.db)

# Now, extract the ENTREZID column from this new data frame
sirt7_entrez_ids <- sirt7_ids_df$ENTREZID

# --- Step 3: Run enrichKEGG with the CORRECTED Entrez ID vectors ---
# Now your inputs are simple vectors of Entrez IDs, and na.omit() will work correctly.

ekegg_7 <- enrichKEGG(gene = na.omit(sirt7_entrez_ids),
                      organism = 'hsa',
                      pvalueCutoff = 1,
                      qvalueCutoff = 0.5,
                      universe = na.omit(background_entrez_ids))


# View the results (if any are found)
head(ekegg_7)

library(enrichplot)
ups_7 <- enrichplot::upsetplot(ekegg_7)
ups_7


# Reactome analysis
BiocManager::install("ReactomePA")
library(ReactomePA)

# The input gene ID should be Entrez gene ID
pathway_results <- enrichPathway(
  gene          = sirt7_entrez_ids,
  organism      = "human",
  universe      = background_entrez_ids,
  pvalueCutoff  = 1,
  qvalueCutoff  = 1,
  pAdjustMethod = "BH",
  readable      = TRUE
)

head(pathway_results)

# View the enriched categories
p7 <- barplot(pathway_results, 
              # showCategory = 7, 
              title = "Reactome analysis for unique SIRT7-KO genes")

# Display the plot
cowplot::plot_grid(p7)







# msigdbr - signature enrichment analysis
# Шаг 1: Анализ с коллекцией Hallmark (H) - Общий обзор
# 1. Загружаем коллекцию Hallmark
hallmark_gsets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

# 2. Запускаем анализ обогащения
enrich_hallmark <- enricher(
  gene = sirt7_entrez_ids,
  universe = background_entrez_ids,
  pvalueCutoff = 1,
  qvalueCutoff = 0.5,
  TERM2GENE = hallmark_gsets
)

# 3. Визуализируем результаты (если они есть)
if (!is.null(enrich_hallmark) && nrow(enrich_hallmark) > 0) {
  p_hallmark <- dotplot(enrich_hallmark, showCategory = 15) +
    ggtitle("Hallmark Pathway Enrichment, SIRT7")
  print(p_hallmark)
} else {
  cat("Не найдено значимых путей в коллекции Hallmark.\n")
}






############### V. GO enrichment analysis using compareCluster function ####

#### 1. For all the DEGs #####


# --- 1. Put all your significant gene data frames into a list ---
# (Assuming sig_1, sig_2, etc. are already loaded and have a 'gene_name' column)
all_sig_tables <- list(
  SIRT1_KO = sig_1,
  SIRT2_KO = sig_2,
  SIRT3_KO = sig_3,
  SIRT4_KO = sig_4,
  SIRT5_KO = sig_5,
  SIRT6_KO = sig_6,
  SIRT7_KO = sig_7
)

# --- 2. Convert the Gene Symbols to Entrez IDs for each list ---
# We use lapply() to apply the same conversion to every table in our list.
# The result will be a named list of Entrez ID vectors.
gene_cluster_entrez <- lapply(all_sig_tables, function(df) {
  # Use bitr for robust ID conversion
  ids_df <- bitr(df$gene_name,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)
  # Return just the vector of Entrez IDs
  return(ids_df$ENTREZID)
})

# Let's check the result: a named list of Entrez IDs
head(gene_cluster_entrez$SIRT1_KO)

# --- 3. Prepare your background (universe) gene list ---
# You correctly specified this should be genes expressed across all experiments.
# The statistically correct way to do this is to take the UNION of all genes
# from your full results tables (res_1_tb, res_2_tb, etc.).
all_res_genes_symbols <- unique(c(
  res_1_tb$gene_name, res_2_tb$gene_name, res_3_tb$gene_name,
  res_4_tb$gene_name, res_5_tb$gene_name, res_6_tb$gene_name,
  res_7_tb$gene_name
))

# Convert the universe symbols to Entrez IDs as well
universe_ids_df <- bitr(na.omit(all_res_genes_symbols),
                        fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)

background_entrez_ids <- universe_ids_df$ENTREZID



# Run the comparative analysis for GO Biological Processes
compare_go_results <- compareCluster(
  geneCluster   = gene_cluster_entrez,
  fun           = "enrichGO",         # Tell it to use the enrichGO function
  universe      = background_entrez_ids,
  OrgDb         = org.Hs.eg.db,       # Argument to pass to enrichGO
  keyType       = "ENTREZID",         # We are providing Entrez IDs
  ont           = "BP",               # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.01
)

# --- Шаг 2: Фильтруем результаты, чтобы оставить только уникальные пути ---

# Проверяем, есть ли вообще результаты для фильтрации
if (!is.null(compare_go_results) && nrow(compare_go_results) > 0) {
  
  # Преобразуем объект в data.frame для удобной работы с dplyr
  results_df <- as.data.frame(compare_go_results)
  
  # Находим уникальные пути. Логика:
  # 1. Группируем по ID онтологии (Description).
  # 2. Считаем, сколько раз встречается каждый ID (n()).
  # 3. Оставляем только те, которые встречаются ровно 1 раз.
  unique_pathways_df <- results_df %>%
    group_by(Description) %>%
    filter(n() == 1) %>%
    ungroup() # Разгруппировываем для дальнейших действий
  
  # --- Шаг 3: Обновляем исходный объект отфильтрованными данными ---
  
  # Создаем копию исходного объекта, чтобы не изменять его
  compare_go_unique <- compare_go_results
  
  # Заменяем внутреннюю таблицу с результатами на нашу отфильтрованную
  compare_go_unique@compareClusterResult <- unique_pathways_df
  
  # --- Шаг 4: Визуализируем только уникальные пути ---
  
  # Проверяем, остались ли пути после фильтрации
  if(nrow(compare_go_unique) > 0) {
    
    # Преобразуем результаты в data.frame для легкого редактирования
    results_to_plot_df <- as.data.frame(compare_go_unique)
    
    # Заменяем длинные имена в колонке 'Cluster' на короткие
    # Функция gsub() отлично подходит для этого: она удаляет "SIRT" и "_KO"
    results_to_plot_df$Cluster <- gsub("SIRT|_KO", "", results_to_plot_df$Cluster)
    
    # Преобразуем колонку Cluster в фактор с правильным порядком уровней
    # Это гарантирует, что на графике они будут идти по порядку от 1 до 7
    results_to_plot_df$Cluster <- factor(results_to_plot_df$Cluster, levels = as.character(1:7))
    
    # Обновляем объект compare_go_unique новыми данными
    compare_go_unique@compareClusterResult <- results_to_plot_df
    
    # --- ТЕПЕРЬ СТРОИМ ГРАФИК ---
    
    p_unique <- dotplot(
      compare_go_unique, 
      showCategory = 2,
      by = "p.adjust"
    ) +
      ggtitle("Unique GO Enrichment Across SIRT KO") +
      xlab("SIRT KO") + # Даем осмысленное имя оси с описаниями
      ylab("")               # И оси с номерами SIRT
    
    print(p_unique)
    
  } else {
    cat("После фильтрации не осталось уникальных путей для отображения.\n")
  }
  
} else {
  cat("В исходных результатах compareCluster нет данных для анализа.\n")
}



# Let's split unique_pathways_df into 7 df for each of the SIRT with its own unique pathways with good qvalue 
list_of_unique_pathway_tables <- split(unique_pathways_df, unique_pathways_df$Cluster)

# Extract the tables with unique pathways 
sirt1_unique_pathways <- list_of_unique_pathway_tables$SIRT1_KO
sirt2_unique_pathways <- list_of_unique_pathway_tables$SIRT2_KO
sirt3_unique_pathways <- list_of_unique_pathway_tables$SIRT3_KO
sirt4_unique_pathways <- list_of_unique_pathway_tables$SIRT4_KO
sirt5_unique_pathways <- list_of_unique_pathway_tables$SIRT5_KO
sirt6_unique_pathways <- list_of_unique_pathway_tables$SIRT6_KO
sirt7_unique_pathways <- list_of_unique_pathway_tables$SIRT7_KO


# 1. SIRT1 
# Check if the data frame is not empty
if (!is.null(sirt1_unique_pathways) && nrow(sirt1_unique_pathways) > 0) {
  
  # Select the top N pathways to plot (e.g., top 10)
  # Arrange by p.adjust to get the most significant ones first
  plot_data <- sirt1_unique_pathways %>%
    arrange(p.adjust) %>%
    slice(1:11) # Take the first 10 rows
  
  # To make ggplot plot them in the correct order, we need to convert
  # the 'Description' column to a factor with the levels in the desired order.
  # We use rev() because ggplot plots from bottom to top.
  plot_data$Description <- factor(plot_data$Description,
                                  levels = rev(plot_data$Description))
  
  # --- Step 4: Create the Bar Plot using ggplot() ---
  
  p_sirt1_barplot_manual <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") + # Use the 'Count' column directly
    scale_fill_continuous(low = "red", high = "blue", name = "p.adjust") +
    labs(
      title = "GO for unique SIRT1_KO genes",
      x = "Count",
      y = NULL # Remove the default "Description" y-axis label
    ) +
    theme_minimal() + # A clean theme to start with
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.text.y = element_text(size = 12)
    )
  
  # --- Step 5: Display the Plot ---
  print(p_sirt1_barplot_manual)
  
} else {
  cat("No unique pathways found for SIRT1-KO to plot.\n")
}


# 2. SIRT2 
# Check if the data frame is not empty
if (!is.null(sirt2_unique_pathways) && nrow(sirt2_unique_pathways) > 0) {
  
  # Select the top N pathways to plot (e.g., top 10)
  # Arrange by p.adjust to get the most significant ones first
  plot_data <- sirt2_unique_pathways %>%
    arrange(p.adjust) %>%
    slice(1:12) # Take the first 10 rows
  
  # To make ggplot plot them in the correct order, we need to convert
  # the 'Description' column to a factor with the levels in the desired order.
  # We use rev() because ggplot plots from bottom to top.
  plot_data$Description <- factor(plot_data$Description,
                                  levels = rev(plot_data$Description))
  
  # --- Step 4: Create the Bar Plot using ggplot() ---
  
  p_sirt2_barplot_manual <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") + # Use the 'Count' column directly
    scale_fill_continuous(low = "red", high = "blue", name = "p.adjust") +
    labs(
      title = "GO for unique SIRT2_KO genes",
      x = "Count",
      y = NULL # Remove the default "Description" y-axis label
    ) +
    theme_minimal() + # A clean theme to start with
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.text.y = element_text(size = 12)
    )
  
  # --- Step 5: Display the Plot ---
  print(p_sirt2_barplot_manual)
  
} else {
  cat("No unique pathways found for SIRT2-KO to plot.\n")
}



# 3. SIRT3 
# Check if the data frame is not empty
if (!is.null(sirt3_unique_pathways) && nrow(sirt3_unique_pathways) > 0) {
  
  # Select the top N pathways to plot (e.g., top 10)
  # Arrange by p.adjust to get the most significant ones first
  plot_data <- sirt3_unique_pathways %>%
    arrange(p.adjust) %>%
    slice(1:10) # Take the first 10 rows
  
  # To make ggplot plot them in the correct order, we need to convert
  # the 'Description' column to a factor with the levels in the desired order.
  # We use rev() because ggplot plots from bottom to top.
  plot_data$Description <- factor(plot_data$Description,
                                  levels = rev(plot_data$Description))
  
  # --- Step 4: Create the Bar Plot using ggplot() ---
  
  p_sirt3_barplot_manual <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") + # Use the 'Count' column directly
    scale_fill_continuous(low = "red", high = "blue", name = "p.adjust") +
    labs(
      title = "GO for unique SIRT3_KO genes",
      x = "Count",
      y = NULL # Remove the default "Description" y-axis label
    ) +
    theme_minimal() + # A clean theme to start with
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.text.y = element_text(size = 12)
    )
  
  # --- Step 5: Display the Plot ---
  print(p_sirt3_barplot_manual)
  
} else {
  cat("No unique pathways found for SIRT3-KO to plot.\n")
}



# 4. SIRT4
# Check if the data frame is not empty
if (!is.null(sirt4_unique_pathways) && nrow(sirt4_unique_pathways) > 0) {
  
  # Select the top N pathways to plot (e.g., top 10)
  # Arrange by p.adjust to get the most significant ones first
  plot_data <- sirt4_unique_pathways %>%
    arrange(p.adjust) %>%
    slice(1:10) # Take the first 10 rows
  
  # To make ggplot plot them in the correct order, we need to convert
  # the 'Description' column to a factor with the levels in the desired order.
  # We use rev() because ggplot plots from bottom to top.
  plot_data$Description <- factor(plot_data$Description,
                                  levels = rev(plot_data$Description))
  
  # --- Step 4: Create the Bar Plot using ggplot() ---
  
  p_sirt4_barplot_manual <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") + # Use the 'Count' column directly
    scale_fill_continuous(low = "red", high = "blue", name = "p.adjust") +
    labs(
      title = "GO for unique SIRT4_KO genes",
      x = "Count",
      y = NULL # Remove the default "Description" y-axis label
    ) +
    theme_minimal() + # A clean theme to start with
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.text.y = element_text(size = 12)
    )
  
  # --- Step 5: Display the Plot ---
  print(p_sirt4_barplot_manual)
  
} else {
  cat("No unique pathways found for SIRT4-KO to plot.\n")
}



# 5. SIRT5
# Check if the data frame is not empty
if (!is.null(sirt5_unique_pathways) && nrow(sirt5_unique_pathways) > 0) {
  
  # Select the top N pathways to plot (e.g., top 10)
  # Arrange by p.adjust to get the most significant ones first
  plot_data <- sirt5_unique_pathways %>%
    arrange(p.adjust) %>%
    slice(1:15) # Take the first 10 rows
  
  # To make ggplot plot them in the correct order, we need to convert
  # the 'Description' column to a factor with the levels in the desired order.
  # We use rev() because ggplot plots from bottom to top.
  plot_data$Description <- factor(plot_data$Description,
                                  levels = rev(plot_data$Description))
  
  # --- Step 4: Create the Bar Plot using ggplot() ---
  
  p_sirt5_barplot_manual <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") + # Use the 'Count' column directly
    scale_fill_continuous(low = "red", high = "blue", name = "p.adjust") +
    labs(
      title = "GO for unique SIRT5_KO genes",
      x = "Count",
      y = NULL # Remove the default "Description" y-axis label
    ) +
    theme_minimal() + # A clean theme to start with
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.text.y = element_text(size = 12)
    )
  
  # --- Step 5: Display the Plot ---
  print(p_sirt5_barplot_manual)
  
} else {
  cat("No unique pathways found for SIRT5-KO to plot.\n")
}



# 6. SIRT6
# Check if the data frame is not empty
if (!is.null(sirt6_unique_pathways) && nrow(sirt6_unique_pathways) > 0) {
  
  # Select the top N pathways to plot (e.g., top 10)
  # Arrange by p.adjust to get the most significant ones first
  plot_data <- sirt6_unique_pathways %>%
    arrange(p.adjust) %>%
    slice(1:10) # Take the first 10 rows
  
  # To make ggplot plot them in the correct order, we need to convert
  # the 'Description' column to a factor with the levels in the desired order.
  # We use rev() because ggplot plots from bottom to top.
  plot_data$Description <- factor(plot_data$Description,
                                  levels = rev(plot_data$Description))
  
  # --- Step 4: Create the Bar Plot using ggplot() ---
  
  p_sirt6_barplot_manual <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") + # Use the 'Count' column directly
    scale_fill_continuous(low = "red", high = "blue", name = "p.adjust") +
    labs(
      title = "GO for unique SIRT6_KO genes",
      x = "Count",
      y = NULL # Remove the default "Description" y-axis label
    ) +
    theme_minimal() + # A clean theme to start with
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.text.y = element_text(size = 12)
    )
  
  # --- Step 5: Display the Plot ---
  print(p_sirt6_barplot_manual)
  
} else {
  cat("No unique pathways found for SIRT6-KO to plot.\n")
}


# 7. SIRT7
# Check if the data frame is not empty
if (!is.null(sirt7_unique_pathways) && nrow(sirt7_unique_pathways) > 0) {
  
  # Select the top N pathways to plot (e.g., top 10)
  # Arrange by p.adjust to get the most significant ones first
  plot_data <- sirt7_unique_pathways %>%
    arrange(p.adjust) %>%
    slice(1:10) # Take the first 10 rows
  
  # To make ggplot plot them in the correct order, we need to convert
  # the 'Description' column to a factor with the levels in the desired order.
  # We use rev() because ggplot plots from bottom to top.
  plot_data$Description <- factor(plot_data$Description,
                                  levels = rev(plot_data$Description))
  
  # --- Step 4: Create the Bar Plot using ggplot() ---
  
  p_sirt7_barplot_manual <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") + # Use the 'Count' column directly
    scale_fill_continuous(low = "red", high = "blue", name = "p.adjust") +
    labs(
      title = "GO for unique SIRT7_KO genes",
      x = "Count",
      y = NULL # Remove the default "Description" y-axis label
    ) +
    theme_minimal() + # A clean theme to start with
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.text.y = element_text(size = 12)
    )
  
  # --- Step 5: Display the Plot ---
  print(p_sirt7_barplot_manual)
  
} else {
  cat("No unique pathways found for SIRT7-KO to plot.\n")
}





# Run KEGG enrichment analysis 
compare_kegg_results <- compareCluster(
  geneCluster = gene_cluster_entrez,
  fun         = "enrichKEGG",         # Tell it to use the enrichKEGG function
  organism    = "hsa",                # Argument to pass to enrichKEGG
  universe    = background_entrez_ids,
  pvalueCutoff = 1,
  qvalueCutoff = 0.5
) # There is no enrichment

# Create a plot of the KEGG comparison results
# This is the standard and most effective visualization
dotplot(
  compare_kegg_results, 
  showCategory = 15,          # Show the top 15 most significant categories
  by = "p.adjust"             # Order them by adjusted p-value
) +
  ggtitle("Comparative KEGG Enrichment Across SIRT KO")




######################## 2. For up-, down-regulated DEGs separately #####

# Create a new list containing only UPREGULATED genes
upregulated_sig_lists <- lapply(all_sig_tables, function(df) {
  
  # Inside the function, we filter the current data frame (df)
  # to keep only rows where log2FoldChange is greater than 0.
  df_filtered <- df %>%
    dplyr::filter(log2FoldChange > 0)
  
  # The function returns the filtered data frame
  return(df_filtered)
})

# Create a new list containing only DOWNREGULATED genes
downregulated_sig_lists <- lapply(all_sig_tables, function(df) {
  
  # Keep only rows where log2FoldChange is less than 0.
  df_filtered <- df %>%
    dplyr::filter(log2FoldChange < 0)
  
  # Return the filtered data frame
  return(df_filtered)
})

# Convert Gene Symbols to Entrez IDs for the UPREGULATED list
upregulated_entrez_lists <- lapply(upregulated_sig_lists, function(df) {
  
  # Check if the data frame is not empty to avoid errors with bitr
  if (nrow(df) > 0) {
    # Use bitr for robust ID conversion from SYMBOL to ENTREZID
    ids_df <- bitr(
      df$gene_name,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    # The function returns only the vector of Entrez IDs
    return(ids_df$ENTREZID)
  } else {
    # If the data frame was empty, return an empty vector
    return(character(0))
  }
})

# Convert Gene Symbols to Entrez IDs for the DOWNREGULATED list
downregulated_entrez_lists <- lapply(downregulated_sig_lists, function(df) {
  
  if (nrow(df) > 0) {
    ids_df <- bitr(
      df$gene_name,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    return(ids_df$ENTREZID)
  } else {
    return(character(0))
  }
})


##### a. For up-regulated genes ###### 

# Run the comparative analysis for GO 
compare_go_results <- compareCluster(
  geneCluster   = upregulated_entrez_lists,
  fun           = "enrichGO",         # Tell it to use the enrichGO function
  universe      = background_entrez_ids,
  OrgDb         = org.Hs.eg.db,       # Argument to pass to enrichGO
  keyType       = "ENTREZID",         # We are providing Entrez IDs
  ont           = "BP",               # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.01
)

# Фильтруем результаты, чтобы оставить только уникальные пути
# Проверяем, есть ли вообще результаты для фильтрации
if (!is.null(compare_go_results) && nrow(compare_go_results) > 0) {
  
  # Преобразуем объект в data.frame для удобной работы с dplyr
  results_df <- as.data.frame(compare_go_results)
  
  # Находим уникальные пути. Логика:
  # 1. Группируем по ID онтологии (Description).
  # 2. Считаем, сколько раз встречается каждый ID (n()).
  # 3. Оставляем только те, которые встречаются ровно 1 раз.
  unique_pathways_df <- results_df %>%
    group_by(Description) %>%
    filter(n() == 1) %>%
    ungroup() # Разгруппировываем для дальнейших действий
  
  # --- Шаг 3: Обновляем исходный объект отфильтрованными данными ---
  
  # Создаем копию исходного объекта, чтобы не изменять его
  compare_go_unique <- compare_go_results
  
  # Заменяем внутреннюю таблицу с результатами на нашу отфильтрованную
  compare_go_unique@compareClusterResult <- unique_pathways_df
  
  # --- Шаг 4: Визуализируем только уникальные пути ---
  
  # Проверяем, остались ли пути после фильтрации
  if(nrow(compare_go_unique) > 0) {
    
    # Преобразуем результаты в data.frame для легкого редактирования
    results_to_plot_df <- as.data.frame(compare_go_unique)
    
    # Заменяем длинные имена в колонке 'Cluster' на короткие
    # Функция gsub() отлично подходит для этого: она удаляет "SIRT" и "_KO"
    results_to_plot_df$Cluster <- gsub("SIRT|_KO", "", results_to_plot_df$Cluster)
    
    # Преобразуем колонку Cluster в фактор с правильным порядком уровней
    # Это гарантирует, что на графике они будут идти по порядку от 1 до 7
    results_to_plot_df$Cluster <- factor(results_to_plot_df$Cluster, levels = as.character(1:7))
    
    # Обновляем объект compare_go_unique новыми данными
    compare_go_unique@compareClusterResult <- results_to_plot_df
    
    # --- ТЕПЕРЬ СТРОИМ ГРАФИК ---
    
    p_unique <- dotplot(
      compare_go_unique, 
      showCategory = 2,
      by = "p.adjust"
    ) +
      ggtitle("Unique GO enrichment across for up- genes") +
      xlab("SIRT KO") + # Даем осмысленное имя оси с описаниями
      ylab("")               # И оси с номерами SIRT
    
    print(p_unique)
    
  } else {
    cat("После фильтрации не осталось уникальных путей для отображения.\n")
  }
  
} else {
  cat("В исходных результатах compareCluster нет данных для анализа.\n")
}

# Split the unique pathways table into a list of tables
list_of_unique_pathway_tables <- split(unique_pathways_df, unique_pathways_df$Cluster)

# Apply simplify() to each SIRT's results and plot



######################## Extract the table for SIRT1
sirt1_unique_df <- list_of_unique_pathway_tables$SIRT1_KO

# Check if there's anything to process
if (!is.null(sirt1_unique_df) && nrow(sirt1_unique_df) > 0) {
  
  # Извлекаем список генов для SIRT1 из исходной таблицы результатов.
  # Это нужно для слота @gene в новом объекте.
  sirt1_gene_list <- gene_cluster_entrez$SIRT1_KO # Предполагается, что у вас есть этот список
  
  # Создаем новый объект класса enrichResult вручную.
  # Мы заполняем его данными из нашей таблицы sirt1_unique_df.
  sirt1_enrich_object_clean <- new("enrichResult",
                                   result = sirt1_unique_df,
                                   pvalueCutoff = 0.05, # Эти значения могут быть любыми
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.1,
                                   gene = as.character(sirt1_gene_list),
                                   universe = as.character(background_entrez_ids),
                                   organism = "HOMO SAPIENS",
                                   keytype = "ENTREZID",
                                   ontology = "BP"
  )
  
  # --- Шаг 4: Упрощаем (simplify) этот новый, чистый объект ---
  sirt1_simplified_results <- simplify(
    sirt1_enrich_object_clean,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
  # --- Шаг 5: Создаем barplot ---
  # Теперь barplot получит объект того класса, который он ожидает.
  
  p_sirt1_barplot <- barplot(
    sirt1_simplified_results,
    showCategory = 15,
    colorBy = "p.adjust",
    title = "GO for unique SIRT1_KO genes"
  ) +
    scale_fill_gradient(low = "red", high = "blue")
  
  # --- Шаг 6: Отображаем график ---
  print(p_sirt1_barplot)
  
} else {
  cat("Не найдено уникальных путей для SIRT1_KO для построения графика.\n")
}

# Преобразуем результаты в data.frame
sirt1_results_df <- as.data.frame(sirt1_simplified_results)

# Проверяем, что в таблице есть строки
if (nrow(sirt1_results_df) > 0) {
  
  # 3. Сохраняем таблицу в CSV файл
  # Этот файл вы будете редактировать вручную
  write.csv(sirt1_results_df, "SIRT1_unique_up_GO_for_grouping.csv", row.names = FALSE)
  
  cat("Таблица 'SIRT1_unique_up_GO_for_grouping.csv' была сохранена в вашей рабочей директории.\n")
  cat("Теперь откройте ее в Excel или Google Sheets для следующего шага.\n")
  
} else {
  cat("Нет данных для анализа и группировки.\n")
}

library(ggplot2)
library(dplyr)
install.packages("tidytext")
library(tidytext)

# --- 1. Читаем вашу новую, сгруппированную таблицу ---
grouped_data <- read.csv2("~/SIRT17_hMSCs/DESeq/SIRT1_unique_up_GO_grouped.csv")

# --- Step 2: Prepare the Data for Plotting ---

# First, ensure the numeric columns are correctly formatted
plot_data <- grouped_data %>%
  mutate(
    # Make sure p.adjust and Count are numeric
    p.adjust = as.numeric(gsub(",", ".", p.adjust)),
    Count = as.numeric(Count),
    # Filter out any rows that failed to parse
    !is.na(p.adjust),
    !is.na(Count)
  )

# For clarity, let's select a limited number of top pathways per theme
plot_data <- plot_data %>%
  group_by(Theme) %>%
  # Arrange by p.adjust within each group and take the top 5
  arrange(p.adjust, .by_group = TRUE) %>%
  slice(1:5) %>%
  ungroup()

# --- CRITICAL STEP for Ordering ---
# To plot correctly, we need to order the 'Description' factor levels
# based on their group (Theme) and their value (Count or p.adjust).
# We also need to order the 'Theme' factor levels.
# Let's order Themes by the max count within each theme for a nice visual.
theme_order <- plot_data %>%
  group_by(Theme) %>%
  summarise(max_count = max(Count)) %>%
  arrange(desc(max_count)) %>%
  pull(Theme)

plot_data <- plot_data %>%
  mutate(
    # Order the main pathway descriptions
    Description = reorder_within(Description, Count, Theme),
    # Order the Theme factor itself
    Theme = factor(Theme, levels = rev(theme_order))
  )


# --- Step 3: Create the Customized Bar Plot ---

p_final <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
  
  # 1. Create the bar plot
  geom_bar(stat = "identity") +
  
  # 2. Add the faceting by Theme. This creates the groups.
  # The 'switch = "y"' argument is the key to moving labels to the right.
  facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  # 3. Apply the correct sorting for the y-axis
  scale_y_reordered() +
  
  # 4. Customize the color scale
  scale_fill_gradientn(
    colors = c("red", "magenta", "blue"), # Your red-to-blue color scheme
    name = "p.adjust"
  ) +
  
  # 5. Customize labels and titles
  labs(
    title = "GO for unique SIRT1_KO genes",
    x = "Count",
    y = NULL # We don't need a y-axis title
  ) +
  
  # 6. Apply a clean theme and customize it to match the reference
  theme_minimal(base_size = 12) +
  theme(
    # Move the group labels (strips) to the outside and align them
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    
    # Remove the background from the group labels
    strip.background = element_blank(),
    
    # Adjust panel spacing
    panel.spacing = unit(0.5, "lines"),
    
    # General plot aesthetics
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines
    panel.grid.minor.x = element_blank()
  )

# --- Step 4: Display the Final Plot ---
print(p_final)



########## Get the number of genes that are found in unique pathways 
# 1. Get the names of the pathways that are unique to SIRT1_KO
sirt1_unique_pathway_names <- unique_pathways_df %>%
  dplyr::filter(Cluster == "SIRT1_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df <- as.data.frame(compare_go_results)

# 3. Filter the full results for SIRT1_KO AND its unique pathways
sirt1_unique_genes_df <- full_results_df %>%
  dplyr::filter(
    Cluster == "SIRT1_KO",
    Description %in% sirt1_unique_pathway_names
  )

# 4. Extract, split, and combine all the gene IDs from the 'geneID' column
# The geneID column is a single string with genes separated by "/"
all_sirt1_unique_genes <- sirt1_unique_genes_df %>%
  pull(geneID) %>%              # Extract the geneID column
  strsplit(split = "/") %>%     # Split each string by the "/"
  unlist() %>%                  # Flatten the list of lists into a single vector
  unique()                      # Keep only the unique gene IDs

# 5. Count them!
number_of_sirt1_genes_in_unique_pathways <- length(all_sirt1_unique_genes)

cat("The number of unique genes contributing to the unique functions of SIRT1 is:",
    number_of_sirt1_genes_in_unique_pathways, "\n")





######################## Extract the table for SIRT2
sirt2_unique_df <- list_of_unique_pathway_tables$SIRT2_KO

# Check if there's anything to process
if (!is.null(sirt2_unique_df) && nrow(sirt2_unique_df) > 0) {
  
  # Извлекаем список генов для SIRT1 из исходной таблицы результатов.
  # Это нужно для слота @gene в новом объекте.
  sirt2_gene_list <- gene_cluster_entrez$SIRT2_KO # Предполагается, что у вас есть этот список
  
  # Создаем новый объект класса enrichResult вручную.
  # Мы заполняем его данными из нашей таблицы sirt1_unique_df.
  sirt2_enrich_object_clean <- new("enrichResult",
                                   result = sirt2_unique_df,
                                   pvalueCutoff = 0.05, # Эти значения могут быть любыми
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.1,
                                   gene = as.character(sirt2_gene_list),
                                   universe = as.character(background_entrez_ids),
                                   organism = "HOMO SAPIENS",
                                   keytype = "ENTREZID",
                                   ontology = "BP"
  )
  
  # --- Шаг 4: Упрощаем (simplify) этот новый, чистый объект ---
  sirt2_simplified_results <- simplify(
    sirt2_enrich_object_clean,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
  # --- Шаг 5: Создаем barplot ---
  # Теперь barplot получит объект того класса, который он ожидает.
  
  p_sirt2_barplot <- barplot(
    sirt2_simplified_results,
    showCategory = 15,
    colorBy = "p.adjust",
    title = "GO for unique SIRT2_KO genes"
  ) +
    scale_fill_gradient(low = "red", high = "blue")
  
  # --- Шаг 6: Отображаем график ---
  print(p_sirt2_barplot)
  
} else {
  cat("Не найдено уникальных путей для SIRT2_KO для построения графика.\n")
}

# Преобразуем результаты в data.frame
sirt2_results_df <- as.data.frame(sirt2_simplified_results)

# Проверяем, что в таблице есть строки
if (nrow(sirt2_results_df) > 0) {
  
  # 3. Сохраняем таблицу в CSV файл
  # Этот файл вы будете редактировать вручную
  write.csv(sirt2_results_df, "SIRT2_unique_up_GO_for_grouping.csv", row.names = FALSE)
  
  cat("Таблица 'SIRT2_unique_up_GO_for_grouping.csv' была сохранена в вашей рабочей директории.\n")
  cat("Теперь откройте ее в Excel или Google Sheets для следующего шага.\n")
  
} else {
  cat("Нет данных для анализа и группировки.\n")
}


# --- 1. Читаем вашу новую, сгруппированную таблицу ---
grouped_data <- read.csv2("~/SIRT17_hMSCs/DESeq/SIRT2_unique_up_GO_grouped.csv")

# --- Step 2: Prepare the Data for Plotting ---

# First, ensure the numeric columns are correctly formatted
plot_data <- grouped_data %>%
  mutate(
    # Make sure p.adjust and Count are numeric
    p.adjust = as.numeric(gsub(",", ".", p.adjust)),
    Count = as.numeric(Count),
    # Filter out any rows that failed to parse
    !is.na(p.adjust),
    !is.na(Count)
  )

# For clarity, let's select a limited number of top pathways per theme
plot_data <- plot_data %>%
  group_by(Theme) %>%
  # Arrange by p.adjust within each group and take the top 5
  arrange(p.adjust, .by_group = TRUE) %>%
  slice(1:5) %>%
  ungroup()

# --- CRITICAL STEP for Ordering ---
# To plot correctly, we need to order the 'Description' factor levels
# based on their group (Theme) and their value (Count or p.adjust).
# We also need to order the 'Theme' factor levels.
# Let's order Themes by the max count within each theme for a nice visual.
theme_order <- plot_data %>%
  group_by(Theme) %>%
  summarise(max_count = max(Count)) %>%
  arrange(desc(max_count)) %>%
  pull(Theme)

plot_data <- plot_data %>%
  mutate(
    # Order the main pathway descriptions
    Description = reorder_within(Description, Count, Theme),
    # Order the Theme factor itself
    Theme = factor(Theme, levels = rev(theme_order))
  )


# --- Step 3: Create the Customized Bar Plot ---

p_final <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
  
  # 1. Create the bar plot
  geom_bar(stat = "identity") +
  
  # 2. Add the faceting by Theme. This creates the groups.
  # The 'switch = "y"' argument is the key to moving labels to the right.
  facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  # 3. Apply the correct sorting for the y-axis
  scale_y_reordered() +
  
  # 4. Customize the color scale
  scale_fill_gradientn(
    colors = c("red", "magenta", "blue"), # Your red-to-blue color scheme
    name = "p.adjust"
  ) +
  
  # 5. Customize labels and titles
  labs(
    title = "GO for unique SIRT2_KO genes",
    x = "Count",
    y = NULL # We don't need a y-axis title
  ) +
  
  # 6. Apply a clean theme and customize it to match the reference
  theme_minimal(base_size = 12) +
  theme(
    # Move the group labels (strips) to the outside and align them
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    
    # Remove the background from the group labels
    strip.background = element_blank(),
    
    # Adjust panel spacing
    panel.spacing = unit(0.5, "lines"),
    
    # General plot aesthetics
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines
    panel.grid.minor.x = element_blank()
  )

# --- Step 4: Display the Final Plot ---
print(p_final)


########## Get the number of genes that are found in unique pathways 
# 1. Get the names of the pathways that are unique to SIRT2_KO
sirt2_unique_pathway_names <- unique_pathways_df %>%
  dplyr::filter(Cluster == "SIRT2_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df <- as.data.frame(compare_go_results)

# 3. Filter the full results for SIRT2_KO AND its unique pathways
sirt2_unique_genes_df <- full_results_df %>%
  dplyr::filter(
    Cluster == "SIRT2_KO",
    Description %in% sirt2_unique_pathway_names
  )

# 4. Extract, split, and combine all the gene IDs from the 'geneID' column
# The geneID column is a single string with genes separated by "/"
all_sirt2_unique_genes <- sirt2_unique_genes_df %>%
  pull(geneID) %>%              # Extract the geneID column
  strsplit(split = "/") %>%     # Split each string by the "/"
  unlist() %>%                  # Flatten the list of lists into a single vector
  unique()                      # Keep only the unique gene IDs

# 5. Count them!
number_of_sirt2_genes_in_unique_pathways <- length(all_sirt2_unique_genes)

cat("The number of unique genes contributing to the unique functions of SIRT2 is:",
    number_of_sirt2_genes_in_unique_pathways, "\n")




######################## Extract the table for SIRT3
sirt3_unique_df <- list_of_unique_pathway_tables$SIRT3_KO

# Check if there's anything to process
if (!is.null(sirt3_unique_df) && nrow(sirt3_unique_df) > 0) {
  
  # Извлекаем список генов для SIRT1 из исходной таблицы результатов.
  # Это нужно для слота @gene в новом объекте.
  sirt3_gene_list <- gene_cluster_entrez$SIRT3_KO # Предполагается, что у вас есть этот список
  
  # Создаем новый объект класса enrichResult вручную.
  # Мы заполняем его данными из нашей таблицы sirt1_unique_df.
  sirt3_enrich_object_clean <- new("enrichResult",
                                   result = sirt3_unique_df,
                                   pvalueCutoff = 0.05, # Эти значения могут быть любыми
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.1,
                                   gene = as.character(sirt3_gene_list),
                                   universe = as.character(background_entrez_ids),
                                   organism = "HOMO SAPIENS",
                                   keytype = "ENTREZID",
                                   ontology = "BP"
  )
  
  # --- Шаг 4: Упрощаем (simplify) этот новый, чистый объект ---
  sirt3_simplified_results <- simplify(
    sirt3_enrich_object_clean,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
  # --- Шаг 5: Создаем barplot ---
  # Теперь barplot получит объект того класса, который он ожидает.
  
  p_sirt3_barplot <- barplot(
    sirt3_simplified_results,
    showCategory = 15,
    colorBy = "p.adjust",
    title = "GO for unique SIRT3_KO genes"
  ) +
    scale_fill_gradient(low = "red", high = "blue")
  
  # --- Шаг 6: Отображаем график ---
  print(p_sirt3_barplot)
  
} else {
  cat("Не найдено уникальных путей для SIRT3_KO для построения графика.\n")
}

# Преобразуем результаты в data.frame
sirt3_results_df <- as.data.frame(sirt3_simplified_results)

# Проверяем, что в таблице есть строки
if (nrow(sirt3_results_df) > 0) {
  
  # 3. Сохраняем таблицу в CSV файл
  # Этот файл вы будете редактировать вручную
  write.csv(sirt3_results_df, "SIRT3_unique_up_GO_for_grouping.csv", row.names = FALSE)
  
  cat("Таблица 'SIRT3_unique_up_GO_for_grouping.csv' была сохранена в вашей рабочей директории.\n")
  cat("Теперь откройте ее в Excel или Google Sheets для следующего шага.\n")
  
} else {
  cat("Нет данных для анализа и группировки.\n")
}


# --- 1. Читаем вашу новую, сгруппированную таблицу ---
grouped_data <- read.csv2("~/SIRT17_hMSCs/DESeq/SIRT3_unique_up_GO_grouped.csv")

# --- Step 2: Prepare the Data for Plotting ---

# First, ensure the numeric columns are correctly formatted
plot_data <- grouped_data %>%
  mutate(
    # Make sure p.adjust and Count are numeric
    p.adjust = as.numeric(gsub(",", ".", p.adjust)),
    Count = as.numeric(Count),
    # Filter out any rows that failed to parse
    !is.na(p.adjust),
    !is.na(Count)
  )

# For clarity, let's select a limited number of top pathways per theme
plot_data <- plot_data %>%
  group_by(Theme) %>%
  # Arrange by p.adjust within each group and take the top 5
  arrange(p.adjust, .by_group = TRUE) %>%
  slice(1:5) %>%
  ungroup()

# --- CRITICAL STEP for Ordering ---
# To plot correctly, we need to order the 'Description' factor levels
# based on their group (Theme) and their value (Count or p.adjust).
# We also need to order the 'Theme' factor levels.
# Let's order Themes by the max count within each theme for a nice visual.
theme_order <- plot_data %>%
  group_by(Theme) %>%
  summarise(max_count = max(Count)) %>%
  arrange(desc(max_count)) %>%
  pull(Theme)

plot_data <- plot_data %>%
  mutate(
    # Order the main pathway descriptions
    Description = reorder_within(Description, Count, Theme),
    # Order the Theme factor itself
    Theme = factor(Theme, levels = rev(theme_order))
  )


# --- Step 3: Create the Customized Bar Plot ---

p_final <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
  
  # 1. Create the bar plot
  geom_bar(stat = "identity") +
  
  # 2. Add the faceting by Theme. This creates the groups.
  # The 'switch = "y"' argument is the key to moving labels to the right.
  facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  # 3. Apply the correct sorting for the y-axis
  scale_y_reordered() +
  
  # 4. Customize the color scale
  scale_fill_gradientn(
    colors = c("red", "magenta", "blue"), # Your red-to-blue color scheme
    name = "p.adjust"
  ) +
  
  # 5. Customize labels and titles
  labs(
    title = "GO for unique SIRT3_KO genes",
    x = "Count",
    y = NULL # We don't need a y-axis title
  ) +
  
  # 6. Apply a clean theme and customize it to match the reference
  theme_minimal(base_size = 12) +
  theme(
    # Move the group labels (strips) to the outside and align them
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    
    # Remove the background from the group labels
    strip.background = element_blank(),
    
    # Adjust panel spacing
    panel.spacing = unit(0.5, "lines"),
    
    # General plot aesthetics
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines
    panel.grid.minor.x = element_blank()
  )

# --- Step 4: Display the Final Plot ---
print(p_final)


########## Get the number of genes that are found in unique pathways 
# 1. Get the names of the pathways that are unique to SIRT3_KO
sirt3_unique_pathway_names <- unique_pathways_df %>%
  dplyr::filter(Cluster == "SIRT3_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df <- as.data.frame(compare_go_results)

# 3. Filter the full results for SIRT3_KO AND its unique pathways
sirt3_unique_genes_df <- full_results_df %>%
  dplyr::filter(
    Cluster == "SIRT3_KO",
    Description %in% sirt3_unique_pathway_names
  )

# 4. Extract, split, and combine all the gene IDs from the 'geneID' column
# The geneID column is a single string with genes separated by "/"
all_sirt3_unique_genes <- sirt3_unique_genes_df %>%
  pull(geneID) %>%              # Extract the geneID column
  strsplit(split = "/") %>%     # Split each string by the "/"
  unlist() %>%                  # Flatten the list of lists into a single vector
  unique()                      # Keep only the unique gene IDs

# 5. Count them!
number_of_sirt3_genes_in_unique_pathways <- length(all_sirt3_unique_genes)

cat("The number of unique genes contributing to the unique functions of SIRT3 is:",
    number_of_sirt3_genes_in_unique_pathways, "\n")





######################## Extract the table for SIRT4
sirt4_unique_df <- list_of_unique_pathway_tables$SIRT4_KO

# Check if there's anything to process
if (!is.null(sirt4_unique_df) && nrow(sirt4_unique_df) > 0) {
  
  # Извлекаем список генов для SIRT1 из исходной таблицы результатов.
  # Это нужно для слота @gene в новом объекте.
  sirt4_gene_list <- gene_cluster_entrez$SIRT4_KO # Предполагается, что у вас есть этот список
  
  # Создаем новый объект класса enrichResult вручную.
  # Мы заполняем его данными из нашей таблицы sirt1_unique_df.
  sirt4_enrich_object_clean <- new("enrichResult",
                                   result = sirt4_unique_df,
                                   pvalueCutoff = 0.05, # Эти значения могут быть любыми
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.1,
                                   gene = as.character(sirt4_gene_list),
                                   universe = as.character(background_entrez_ids),
                                   organism = "HOMO SAPIENS",
                                   keytype = "ENTREZID",
                                   ontology = "BP"
  )
  
  # --- Шаг 4: Упрощаем (simplify) этот новый, чистый объект ---
  sirt4_simplified_results <- simplify(
    sirt4_enrich_object_clean,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
  # --- Шаг 5: Создаем barplot ---
  # Теперь barplot получит объект того класса, который он ожидает.
  
  p_sirt4_barplot <- barplot(
    sirt4_simplified_results,
    showCategory = 15,
    colorBy = "p.adjust",
    title = "GO for unique SIRT4_KO genes"
  ) +
    scale_fill_gradient(low = "red", high = "blue")
  
  # --- Шаг 6: Отображаем график ---
  print(p_sirt4_barplot)
  
} else {
  cat("Не найдено уникальных путей для SIRT4_KO для построения графика.\n")
}

# Преобразуем результаты в data.frame
sirt4_results_df <- as.data.frame(sirt4_simplified_results)

# Проверяем, что в таблице есть строки
if (nrow(sirt4_results_df) > 0) {
  
  # 3. Сохраняем таблицу в CSV файл
  # Этот файл вы будете редактировать вручную
  write.csv(sirt4_results_df, "SIRT4_unique_up_GO_for_grouping.csv", row.names = FALSE)
  
  cat("Таблица 'SIRT4_unique_up_GO_for_grouping.csv' была сохранена в вашей рабочей директории.\n")
  cat("Теперь откройте ее в Excel или Google Sheets для следующего шага.\n")
  
} else {
  cat("Нет данных для анализа и группировки.\n")
}


# --- 1. Читаем вашу новую, сгруппированную таблицу ---
grouped_data <- read.csv2("~/SIRT17_hMSCs/DESeq/SIRT4_unique_up_GO_grouped.csv")

# --- Step 2: Prepare the Data for Plotting ---

# First, ensure the numeric columns are correctly formatted
plot_data <- grouped_data %>%
  mutate(
    # Make sure p.adjust and Count are numeric
    p.adjust = as.numeric(gsub(",", ".", p.adjust)),
    Count = as.numeric(Count),
    # Filter out any rows that failed to parse
    !is.na(p.adjust),
    !is.na(Count)
  )

# For clarity, let's select a limited number of top pathways per theme
plot_data <- plot_data %>%
  group_by(Theme) %>%
  # Arrange by p.adjust within each group and take the top 5
  arrange(p.adjust, .by_group = TRUE) %>%
  slice(1:5) %>%
  ungroup()

# --- CRITICAL STEP for Ordering ---
# To plot correctly, we need to order the 'Description' factor levels
# based on their group (Theme) and their value (Count or p.adjust).
# We also need to order the 'Theme' factor levels.
# Let's order Themes by the max count within each theme for a nice visual.
theme_order <- plot_data %>%
  group_by(Theme) %>%
  summarise(max_count = max(Count)) %>%
  arrange(desc(max_count)) %>%
  pull(Theme)

plot_data <- plot_data %>%
  mutate(
    # Order the main pathway descriptions
    Description = reorder_within(Description, Count, Theme),
    # Order the Theme factor itself
    Theme = factor(Theme, levels = rev(theme_order))
  )


# --- Step 3: Create the Customized Bar Plot ---

p_final <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
  
  # 1. Create the bar plot
  geom_bar(stat = "identity") +
  
  # 2. Add the faceting by Theme. This creates the groups.
  # The 'switch = "y"' argument is the key to moving labels to the right.
  facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  # 3. Apply the correct sorting for the y-axis
  scale_y_reordered() +
  
  # 4. Customize the color scale
  scale_fill_gradientn(
    colors = c("red", "magenta", "blue"), # Your red-to-blue color scheme
    name = "p.adjust"
  ) +
  
  # 5. Customize labels and titles
  labs(
    title = "GO for unique SIRT4_KO genes",
    x = "Count",
    y = NULL # We don't need a y-axis title
  ) +
  
  # 6. Apply a clean theme and customize it to match the reference
  theme_minimal(base_size = 12) +
  theme(
    # Move the group labels (strips) to the outside and align them
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    
    # Remove the background from the group labels
    strip.background = element_blank(),
    
    # Adjust panel spacing
    panel.spacing = unit(0.5, "lines"),
    
    # General plot aesthetics
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines
    panel.grid.minor.x = element_blank()
  )

# --- Step 4: Display the Final Plot ---
print(p_final)



########## Get the number of genes that are found in unique pathways 
# 1. Get the names of the pathways that are unique to SIRT4_KO
sirt4_unique_pathway_names <- unique_pathways_df %>%
  dplyr::filter(Cluster == "SIRT4_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df <- as.data.frame(compare_go_results)

# 3. Filter the full results for SIRT4_KO AND its unique pathways
sirt4_unique_genes_df <- full_results_df %>%
  dplyr::filter(
    Cluster == "SIRT4_KO",
    Description %in% sirt4_unique_pathway_names
  )

# 4. Extract, split, and combine all the gene IDs from the 'geneID' column
# The geneID column is a single string with genes separated by "/"
all_sirt4_unique_genes <- sirt4_unique_genes_df %>%
  pull(geneID) %>%              # Extract the geneID column
  strsplit(split = "/") %>%     # Split each string by the "/"
  unlist() %>%                  # Flatten the list of lists into a single vector
  unique()                      # Keep only the unique gene IDs

# 5. Count them!
number_of_sirt4_genes_in_unique_pathways <- length(all_sirt4_unique_genes)

cat("The number of unique genes contributing to the unique functions of SIRT4 is:",
    number_of_sirt4_genes_in_unique_pathways, "\n")




######################## Extract the table for SIRT5
sirt5_unique_df <- list_of_unique_pathway_tables$SIRT5_KO

# Check if there's anything to process
if (!is.null(sirt5_unique_df) && nrow(sirt5_unique_df) > 0) {
  
  # Извлекаем список генов для SIRT1 из исходной таблицы результатов.
  # Это нужно для слота @gene в новом объекте.
  sirt5_gene_list <- gene_cluster_entrez$SIRT5_KO # Предполагается, что у вас есть этот список
  
  # Создаем новый объект класса enrichResult вручную.
  # Мы заполняем его данными из нашей таблицы sirt1_unique_df.
  sirt5_enrich_object_clean <- new("enrichResult",
                                   result = sirt5_unique_df,
                                   pvalueCutoff = 0.05, # Эти значения могут быть любыми
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.1,
                                   gene = as.character(sirt5_gene_list),
                                   universe = as.character(background_entrez_ids),
                                   organism = "HOMO SAPIENS",
                                   keytype = "ENTREZID",
                                   ontology = "BP"
  )
  
  # --- Шаг 4: Упрощаем (simplify) этот новый, чистый объект ---
  sirt5_simplified_results <- simplify(
    sirt5_enrich_object_clean,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
  # --- Шаг 5: Создаем barplot ---
  # Теперь barplot получит объект того класса, который он ожидает.
  
  p_sirt5_barplot <- barplot(
    sirt5_simplified_results,
    showCategory = 15,
    colorBy = "p.adjust",
    title = "GO for unique SIRT5_KO genes"
  ) +
    scale_fill_gradient(low = "red", high = "blue")
  
  # --- Шаг 6: Отображаем график ---
  print(p_sirt5_barplot)
  
} else {
  cat("Не найдено уникальных путей для SIRT5_KO для построения графика.\n")
}

# Преобразуем результаты в data.frame
sirt5_results_df <- as.data.frame(sirt5_simplified_results)

# Проверяем, что в таблице есть строки
if (nrow(sirt5_results_df) > 0) {
  
  # 3. Сохраняем таблицу в CSV файл
  # Этот файл вы будете редактировать вручную
  write.csv(sirt5_results_df, "SIRT5_unique_up_GO_for_grouping.csv", row.names = FALSE)
  
  cat("Таблица 'SIRT5_unique_up_GO_for_grouping.csv' была сохранена в вашей рабочей директории.\n")
  cat("Теперь откройте ее в Excel или Google Sheets для следующего шага.\n")
  
} else {
  cat("Нет данных для анализа и группировки.\n")
}


# --- 1. Читаем вашу новую, сгруппированную таблицу ---
grouped_data <- read.csv2("~/SIRT17_hMSCs/DESeq/SIRT5_unique_up_GO_grouped.csv")

# --- Step 2: Prepare the Data for Plotting ---

# First, ensure the numeric columns are correctly formatted
plot_data <- grouped_data %>%
  mutate(
    # Make sure p.adjust and Count are numeric
    p.adjust = as.numeric(gsub(",", ".", p.adjust)),
    Count = as.numeric(Count),
    # Filter out any rows that failed to parse
    !is.na(p.adjust),
    !is.na(Count)
  )

# For clarity, let's select a limited number of top pathways per theme
plot_data <- plot_data %>%
  group_by(Theme) %>%
  # Arrange by p.adjust within each group and take the top 5
  arrange(p.adjust, .by_group = TRUE) %>%
  slice(1:5) %>%
  ungroup()

# --- CRITICAL STEP for Ordering ---
# To plot correctly, we need to order the 'Description' factor levels
# based on their group (Theme) and their value (Count or p.adjust).
# We also need to order the 'Theme' factor levels.
# Let's order Themes by the max count within each theme for a nice visual.
theme_order <- plot_data %>%
  group_by(Theme) %>%
  summarise(max_count = max(Count)) %>%
  arrange(desc(max_count)) %>%
  pull(Theme)

plot_data <- plot_data %>%
  mutate(
    # Order the main pathway descriptions
    Description = reorder_within(Description, Count, Theme),
    # Order the Theme factor itself
    Theme = factor(Theme, levels = rev(theme_order))
  )


# --- Step 3: Create the Customized Bar Plot ---

p_final <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
  
  # 1. Create the bar plot
  geom_bar(stat = "identity") +
  
  # 2. Add the faceting by Theme. This creates the groups.
  # The 'switch = "y"' argument is the key to moving labels to the right.
  facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  # 3. Apply the correct sorting for the y-axis
  scale_y_reordered() +
  
  # 4. Customize the color scale
  scale_fill_gradientn(
    colors = c("red", "magenta", "blue"), # Your red-to-blue color scheme
    name = "p.adjust"
  ) +
  
  # 5. Customize labels and titles
  labs(
    title = "GO for unique SIRT5_KO genes",
    x = "Count",
    y = NULL # We don't need a y-axis title
  ) +
  
  # 6. Apply a clean theme and customize it to match the reference
  theme_minimal(base_size = 12) +
  theme(
    # Move the group labels (strips) to the outside and align them
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    
    # Remove the background from the group labels
    strip.background = element_blank(),
    
    # Adjust panel spacing
    panel.spacing = unit(0.5, "lines"),
    
    # General plot aesthetics
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines
    panel.grid.minor.x = element_blank()
  )

# --- Step 4: Display the Final Plot ---
print(p_final)



########## Get the number of genes that are found in unique pathways 
# 1. Get the names of the pathways that are unique to SIRT5_KO
sirt5_unique_pathway_names <- unique_pathways_df %>%
  dplyr::filter(Cluster == "SIRT5_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df <- as.data.frame(compare_go_results)

# 3. Filter the full results for SIRT5_KO AND its unique pathways
sirt5_unique_genes_df <- full_results_df %>%
  dplyr::filter(
    Cluster == "SIRT5_KO",
    Description %in% sirt5_unique_pathway_names
  )

# 4. Extract, split, and combine all the gene IDs from the 'geneID' column
# The geneID column is a single string with genes separated by "/"
all_sirt5_unique_genes <- sirt5_unique_genes_df %>%
  pull(geneID) %>%              # Extract the geneID column
  strsplit(split = "/") %>%     # Split each string by the "/"
  unlist() %>%                  # Flatten the list of lists into a single vector
  unique()                      # Keep only the unique gene IDs

# 5. Count them!
number_of_sirt5_genes_in_unique_pathways <- length(all_sirt5_unique_genes)

cat("The number of unique genes contributing to the unique functions of SIRT5 is:",
    number_of_sirt5_genes_in_unique_pathways, "\n")



######################## Extract the table for SIRT6
sirt6_unique_df <- list_of_unique_pathway_tables$SIRT6_KO

# Check if there's anything to process
if (!is.null(sirt6_unique_df) && nrow(sirt6_unique_df) > 0) {
  
  # Извлекаем список генов для SIRT1 из исходной таблицы результатов.
  # Это нужно для слота @gene в новом объекте.
  sirt6_gene_list <- gene_cluster_entrez$SIRT6_KO # Предполагается, что у вас есть этот список
  
  # Создаем новый объект класса enrichResult вручную.
  # Мы заполняем его данными из нашей таблицы sirt1_unique_df.
  sirt6_enrich_object_clean <- new("enrichResult",
                                   result = sirt6_unique_df,
                                   pvalueCutoff = 0.05, # Эти значения могут быть любыми
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.1,
                                   gene = as.character(sirt6_gene_list),
                                   universe = as.character(background_entrez_ids),
                                   organism = "HOMO SAPIENS",
                                   keytype = "ENTREZID",
                                   ontology = "BP"
  )
  
  # --- Шаг 4: Упрощаем (simplify) этот новый, чистый объект ---
  sirt6_simplified_results <- simplify(
    sirt6_enrich_object_clean,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
  # --- Шаг 5: Создаем barplot ---
  # Теперь barplot получит объект того класса, который он ожидает.
  
  p_sirt6_barplot <- barplot(
    sirt6_simplified_results,
    showCategory = 15,
    colorBy = "p.adjust",
    title = "GO for unique SIRT6_KO genes"
  ) +
    scale_fill_gradient(low = "red", high = "blue")
  
  # --- Шаг 6: Отображаем график ---
  print(p_sirt6_barplot)
  
} else {
  cat("Не найдено уникальных путей для SIRT6_KO для построения графика.\n")
}

# Преобразуем результаты в data.frame
sirt6_results_df <- as.data.frame(sirt6_simplified_results)

# Проверяем, что в таблице есть строки
if (nrow(sirt6_results_df) > 0) {
  
  # 3. Сохраняем таблицу в CSV файл
  # Этот файл вы будете редактировать вручную
  write.csv(sirt6_results_df, "SIRT6_unique_up_GO_for_grouping.csv", row.names = FALSE)
  
  cat("Таблица 'SIRT6_unique_up_GO_for_grouping.csv' была сохранена в вашей рабочей директории.\n")
  cat("Теперь откройте ее в Excel или Google Sheets для следующего шага.\n")
  
} else {
  cat("Нет данных для анализа и группировки.\n")
}


# --- 1. Читаем вашу новую, сгруппированную таблицу ---
grouped_data <- read.csv2("~/SIRT17_hMSCs/DESeq/SIRT6_unique_up_GO_grouped.csv")

# --- Step 2: Prepare the Data for Plotting ---

# First, ensure the numeric columns are correctly formatted
plot_data <- grouped_data %>%
  mutate(
    # Make sure p.adjust and Count are numeric
    p.adjust = as.numeric(gsub(",", ".", p.adjust)),
    Count = as.numeric(Count),
    # Filter out any rows that failed to parse
    !is.na(p.adjust),
    !is.na(Count)
  )

# For clarity, let's select a limited number of top pathways per theme
plot_data <- plot_data %>%
  group_by(Theme) %>%
  # Arrange by p.adjust within each group and take the top 5
  arrange(p.adjust, .by_group = TRUE) %>%
  slice(1:5) %>%
  ungroup()

# --- CRITICAL STEP for Ordering ---
# To plot correctly, we need to order the 'Description' factor levels
# based on their group (Theme) and their value (Count or p.adjust).
# We also need to order the 'Theme' factor levels.
# Let's order Themes by the max count within each theme for a nice visual.
theme_order <- plot_data %>%
  group_by(Theme) %>%
  summarise(max_count = max(Count)) %>%
  arrange(desc(max_count)) %>%
  pull(Theme)

plot_data <- plot_data %>%
  mutate(
    # Order the main pathway descriptions
    Description = reorder_within(Description, Count, Theme),
    # Order the Theme factor itself
    Theme = factor(Theme, levels = rev(theme_order))
  )


# --- Step 3: Create the Customized Bar Plot ---

p_final <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
  
  # 1. Create the bar plot
  geom_bar(stat = "identity") +
  
  # 2. Add the faceting by Theme. This creates the groups.
  # The 'switch = "y"' argument is the key to moving labels to the right.
  facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  # 3. Apply the correct sorting for the y-axis
  scale_y_reordered() +
  
  # 4. Customize the color scale
  scale_fill_gradientn(
    colors = c("red", "magenta", "blue"), # Your red-to-blue color scheme
    name = "p.adjust"
  ) +
  
  # 5. Customize labels and titles
  labs(
    title = "GO for unique SIRT6_KO genes",
    x = "Count",
    y = NULL # We don't need a y-axis title
  ) +
  
  # 6. Apply a clean theme and customize it to match the reference
  theme_minimal(base_size = 12) +
  theme(
    # Move the group labels (strips) to the outside and align them
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    
    # Remove the background from the group labels
    strip.background = element_blank(),
    
    # Adjust panel spacing
    panel.spacing = unit(0.5, "lines"),
    
    # General plot aesthetics
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines
    panel.grid.minor.x = element_blank()
  )

# --- Step 4: Display the Final Plot ---
print(p_final)



########## Get the number of genes that are found in unique pathways 
# 1. Get the names of the pathways that are unique to SIRT6_KO
sirt6_unique_pathway_names <- unique_pathways_df %>%
  dplyr::filter(Cluster == "SIRT6_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df <- as.data.frame(compare_go_results)

# 3. Filter the full results for SIRT6_KO AND its unique pathways
sirt6_unique_genes_df <- full_results_df %>%
  dplyr::filter(
    Cluster == "SIRT6_KO",
    Description %in% sirt6_unique_pathway_names
  )

# 4. Extract, split, and combine all the gene IDs from the 'geneID' column
# The geneID column is a single string with genes separated by "/"
all_sirt6_unique_genes <- sirt6_unique_genes_df %>%
  pull(geneID) %>%              # Extract the geneID column
  strsplit(split = "/") %>%     # Split each string by the "/"
  unlist() %>%                  # Flatten the list of lists into a single vector
  unique()                      # Keep only the unique gene IDs

# 5. Count them!
number_of_sirt6_genes_in_unique_pathways <- length(all_sirt6_unique_genes)

cat("The number of unique genes contributing to the unique functions of SIRT6 is:",
    number_of_sirt6_genes_in_unique_pathways, "\n")




######################## Extract the table for SIRT7
sirt7_unique_df <- list_of_unique_pathway_tables$SIRT7_KO

# Check if there's anything to process
if (!is.null(sirt7_unique_df) && nrow(sirt7_unique_df) > 0) {
  
  # Извлекаем список генов для SIRT1 из исходной таблицы результатов.
  # Это нужно для слота @gene в новом объекте.
  sirt7_gene_list <- gene_cluster_entrez$SIRT7_KO # Предполагается, что у вас есть этот список
  
  # Создаем новый объект класса enrichResult вручную.
  # Мы заполняем его данными из нашей таблицы sirt1_unique_df.
  sirt7_enrich_object_clean <- new("enrichResult",
                                   result = sirt7_unique_df,
                                   pvalueCutoff = 0.05, # Эти значения могут быть любыми
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.1,
                                   gene = as.character(sirt7_gene_list),
                                   universe = as.character(background_entrez_ids),
                                   organism = "HOMO SAPIENS",
                                   keytype = "ENTREZID",
                                   ontology = "BP"
  )
  
  # --- Шаг 4: Упрощаем (simplify) этот новый, чистый объект ---
  sirt7_simplified_results <- simplify(
    sirt7_enrich_object_clean,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
  # --- Шаг 5: Создаем barplot ---
  # Теперь barplot получит объект того класса, который он ожидает.
  
  p_sirt7_barplot <- barplot(
    sirt7_simplified_results,
    showCategory = 15,
    colorBy = "p.adjust",
    title = "GO for unique SIRT7_KO genes"
  ) +
    scale_fill_gradient(low = "red", high = "blue")
  
  # --- Шаг 6: Отображаем график ---
  print(p_sirt7_barplot)
  
} else {
  cat("Не найдено уникальных путей для SIRT7_KO для построения графика.\n")
}

# Преобразуем результаты в data.frame
sirt7_results_df <- as.data.frame(sirt7_simplified_results)

# Проверяем, что в таблице есть строки
if (nrow(sirt7_results_df) > 0) {
  
  # 3. Сохраняем таблицу в CSV файл
  # Этот файл вы будете редактировать вручную
  write.csv(sirt7_results_df, "SIRT7_unique_up_GO_for_grouping.csv", row.names = FALSE)
  
  cat("Таблица 'SIRT7_unique_up_GO_for_grouping.csv' была сохранена в вашей рабочей директории.\n")
  cat("Теперь откройте ее в Excel или Google Sheets для следующего шага.\n")
  
} else {
  cat("Нет данных для анализа и группировки.\n")
}


# --- 1. Читаем вашу новую, сгруппированную таблицу ---
grouped_data <- read.csv2("~/SIRT17_hMSCs/DESeq/SIRT7_unique_up_GO_grouped.csv")

# --- Step 2: Prepare the Data for Plotting ---

# First, ensure the numeric columns are correctly formatted
plot_data <- grouped_data %>%
  mutate(
    # Make sure p.adjust and Count are numeric
    p.adjust = as.numeric(gsub(",", ".", p.adjust)),
    Count = as.numeric(Count),
    # Filter out any rows that failed to parse
    !is.na(p.adjust),
    !is.na(Count)
  )

# For clarity, let's select a limited number of top pathways per theme
plot_data <- plot_data %>%
  group_by(Theme) %>%
  # Arrange by p.adjust within each group and take the top 5
  arrange(p.adjust, .by_group = TRUE) %>%
  slice(1:5) %>%
  ungroup()

# --- CRITICAL STEP for Ordering ---
# To plot correctly, we need to order the 'Description' factor levels
# based on their group (Theme) and their value (Count or p.adjust).
# We also need to order the 'Theme' factor levels.
# Let's order Themes by the max count within each theme for a nice visual.
theme_order <- plot_data %>%
  group_by(Theme) %>%
  summarise(max_count = max(Count)) %>%
  arrange(desc(max_count)) %>%
  pull(Theme)

plot_data <- plot_data %>%
  mutate(
    # Order the main pathway descriptions
    Description = reorder_within(Description, Count, Theme),
    # Order the Theme factor itself
    Theme = factor(Theme, levels = rev(theme_order))
  )


# --- Step 3: Create the Customized Bar Plot ---

p_final <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
  
  # 1. Create the bar plot
  geom_bar(stat = "identity") +
  
  # 2. Add the faceting by Theme. This creates the groups.
  # The 'switch = "y"' argument is the key to moving labels to the right.
  facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  # 3. Apply the correct sorting for the y-axis
  scale_y_reordered() +
  
  # 4. Customize the color scale
  scale_fill_gradientn(
    colors = c("red", "magenta", "blue"), # Your red-to-blue color scheme
    name = "p.adjust"
  ) +
  
  # 5. Customize labels and titles
  labs(
    title = "GO for unique SIRT6_KO genes",
    x = "Count",
    y = NULL # We don't need a y-axis title
  ) +
  
  # 6. Apply a clean theme and customize it to match the reference
  theme_minimal(base_size = 12) +
  theme(
    # Move the group labels (strips) to the outside and align them
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    
    # Remove the background from the group labels
    strip.background = element_blank(),
    
    # Adjust panel spacing
    panel.spacing = unit(0.5, "lines"),
    
    # General plot aesthetics
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines
    panel.grid.minor.x = element_blank()
  )

# --- Step 4: Display the Final Plot ---
print(p_final)



########## Get the number of genes that are found in unique pathways 
# 1. Get the names of the pathways that are unique to SIRT7_KO
sirt7_unique_pathway_names <- unique_pathways_df %>%
  dplyr::filter(Cluster == "SIRT7_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df <- as.data.frame(compare_go_results)

# 3. Filter the full results for SIRT7_KO AND its unique pathways
sirt7_unique_genes_df <- full_results_df %>%
  dplyr::filter(
    Cluster == "SIRT7_KO",
    Description %in% sirt7_unique_pathway_names
  )

# 4. Extract, split, and combine all the gene IDs from the 'geneID' column
# The geneID column is a single string with genes separated by "/"
all_sirt7_unique_genes <- sirt7_unique_genes_df %>%
  pull(geneID) %>%              # Extract the geneID column
  strsplit(split = "/") %>%     # Split each string by the "/"
  unlist() %>%                  # Flatten the list of lists into a single vector
  unique()                      # Keep only the unique gene IDs

# 5. Count them!
number_of_sirt7_genes_in_unique_pathways <- length(all_sirt7_unique_genes)

cat("The number of unique genes contributing to the unique functions of SIRT7 is:",
    number_of_sirt7_genes_in_unique_pathways, "\n")






####################### b. For down-regulated genes ######

# Run the comparative analysis for GO 
compare_go_results <- compareCluster(
  geneCluster   = downregulated_entrez_lists,
  fun           = "enrichGO",         # Tell it to use the enrichGO function
  universe      = background_entrez_ids,
  OrgDb         = org.Hs.eg.db,       # Argument to pass to enrichGO
  keyType       = "ENTREZID",         # We are providing Entrez IDs
  ont           = "BP",               # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.01
)

# Фильтруем результаты, чтобы оставить только уникальные пути 
# Проверяем, есть ли вообще результаты для фильтрации
if (!is.null(compare_go_results) && nrow(compare_go_results) > 0) {
  
  # Преобразуем объект в data.frame для удобной работы с dplyr
  results_df <- as.data.frame(compare_go_results)
  
  # Находим уникальные пути. Логика:
  # 1. Группируем по ID онтологии (Description).
  # 2. Считаем, сколько раз встречается каждый ID (n()).
  # 3. Оставляем только те, которые встречаются ровно 1 раз.
  unique_pathways_df <- results_df %>%
    group_by(Description) %>%
    filter(n() == 1) %>%
    ungroup() # Разгруппировываем для дальнейших действий
  
  # --- Шаг 3: Обновляем исходный объект отфильтрованными данными ---
  
  # Создаем копию исходного объекта, чтобы не изменять его
  compare_go_unique <- compare_go_results
  
  # Заменяем внутреннюю таблицу с результатами на нашу отфильтрованную
  compare_go_unique@compareClusterResult <- unique_pathways_df
  
  # --- Шаг 4: Визуализируем только уникальные пути ---
  
  # Проверяем, остались ли пути после фильтрации
  if(nrow(compare_go_unique) > 0) {
    
    # Преобразуем результаты в data.frame для легкого редактирования
    results_to_plot_df <- as.data.frame(compare_go_unique)
    
    # Заменяем длинные имена в колонке 'Cluster' на короткие
    # Функция gsub() отлично подходит для этого: она удаляет "SIRT" и "_KO"
    results_to_plot_df$Cluster <- gsub("SIRT|_KO", "", results_to_plot_df$Cluster)
    
    # Преобразуем колонку Cluster в фактор с правильным порядком уровней
    # Это гарантирует, что на графике они будут идти по порядку от 1 до 7
    results_to_plot_df$Cluster <- factor(results_to_plot_df$Cluster, levels = as.character(1:7))
    
    # Обновляем объект compare_go_unique новыми данными
    compare_go_unique@compareClusterResult <- results_to_plot_df
    
    # --- ТЕПЕРЬ СТРОИМ ГРАФИК ---
    
    p_unique <- dotplot(
      compare_go_unique, 
      showCategory = 2,
      by = "p.adjust"
    ) +
      ggtitle("Unique GO enrichment across for down- genes") +
      xlab("SIRT KO") + # Даем осмысленное имя оси с описаниями
      ylab("")               # И оси с номерами SIRT
    
    print(p_unique)
    
  } else {
    cat("После фильтрации не осталось уникальных путей для отображения.\n")
  }
  
} else {
  cat("В исходных результатах compareCluster нет данных для анализа.\n")
}


# Split the unique pathways table into a list of tables
list_of_unique_pathway_tables <- split(unique_pathways_df, unique_pathways_df$Cluster)

# Apply simplify() to each SIRT's results and plot


######################## Extract the table for SIRT1
sirt1_unique_df <- list_of_unique_pathway_tables$SIRT1_KO

# Check if there's anything to process
if (!is.null(sirt1_unique_df) && nrow(sirt1_unique_df) > 0) {
  
  # Извлекаем список генов для SIRT1 из исходной таблицы результатов.
  # Это нужно для слота @gene в новом объекте.
  sirt1_gene_list <- gene_cluster_entrez$SIRT1_KO # Предполагается, что у вас есть этот список
  
  # Создаем новый объект класса enrichResult вручную.
  # Мы заполняем его данными из нашей таблицы sirt1_unique_df.
  sirt1_enrich_object_clean <- new("enrichResult",
                                   result = sirt1_unique_df,
                                   pvalueCutoff = 0.05, # Эти значения могут быть любыми
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.1,
                                   gene = as.character(sirt1_gene_list),
                                   universe = as.character(background_entrez_ids),
                                   organism = "HOMO SAPIENS",
                                   keytype = "ENTREZID",
                                   ontology = "BP"
  )
  
  # --- Шаг 4: Упрощаем (simplify) этот новый, чистый объект ---
  sirt1_simplified_results <- simplify(
    sirt1_enrich_object_clean,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
  # --- Шаг 5: Создаем barplot ---
  # Теперь barplot получит объект того класса, который он ожидает.
  
  p_sirt1_barplot <- barplot(
    sirt1_simplified_results,
    showCategory = 15,
    colorBy = "p.adjust",
    title = "GO for unique SIRT1_KO genes"
  ) +
    scale_fill_gradient(low = "red", high = "blue")
  
  # --- Шаг 6: Отображаем график ---
  print(p_sirt1_barplot)
  
} else {
  cat("Не найдено уникальных путей для SIRT1_KO для построения графика.\n")
}

# Преобразуем результаты в data.frame
sirt1_results_df <- as.data.frame(sirt1_simplified_results)

# Проверяем, что в таблице есть строки
if (nrow(sirt1_results_df) > 0) {
  
  # 3. Сохраняем таблицу в CSV файл
  # Этот файл вы будете редактировать вручную
  write.csv(sirt1_results_df, "SIRT1_unique_down_GO_for_grouping.csv", row.names = FALSE)
  
  cat("Таблица 'SIRT1_unique_down_GO_for_grouping.csv' была сохранена в вашей рабочей директории.\n")
  cat("Теперь откройте ее в Excel или Google Sheets для следующего шага.\n")
  
} else {
  cat("Нет данных для анализа и группировки.\n")
}


# --- 1. Читаем вашу новую, сгруппированную таблицу ---
grouped_data <- read.csv2("~/SIRT17_hMSCs/DESeq/SIRT1_unique_down_GO_grouped.csv")

# --- Step 2: Prepare the Data for Plotting ---

# First, ensure the numeric columns are correctly formatted
plot_data <- grouped_data %>%
  mutate(
    # Make sure p.adjust and Count are numeric
    p.adjust = as.numeric(gsub(",", ".", p.adjust)),
    Count = as.numeric(Count),
    # Filter out any rows that failed to parse
    !is.na(p.adjust),
    !is.na(Count)
  )

# For clarity, let's select a limited number of top pathways per theme
plot_data <- plot_data %>%
  group_by(Theme) %>%
  # Arrange by p.adjust within each group and take the top 5
  arrange(p.adjust, .by_group = TRUE) %>%
  slice(1:5) %>%
  ungroup()

# --- CRITICAL STEP for Ordering ---
# To plot correctly, we need to order the 'Description' factor levels
# based on their group (Theme) and their value (Count or p.adjust).
# We also need to order the 'Theme' factor levels.
# Let's order Themes by the max count within each theme for a nice visual.
theme_order <- plot_data %>%
  group_by(Theme) %>%
  summarise(max_count = max(Count)) %>%
  arrange(desc(max_count)) %>%
  pull(Theme)

plot_data <- plot_data %>%
  mutate(
    # Order the main pathway descriptions
    Description = reorder_within(Description, Count, Theme),
    # Order the Theme factor itself
    Theme = factor(Theme, levels = rev(theme_order))
  )


# --- Step 3: Create the Customized Bar Plot ---

p_final <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
  
  # 1. Create the bar plot
  geom_bar(stat = "identity") +
  
  # 2. Add the faceting by Theme. This creates the groups.
  # The 'switch = "y"' argument is the key to moving labels to the right.
  facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  # 3. Apply the correct sorting for the y-axis
  scale_y_reordered() +
  
  # 4. Customize the color scale
  scale_fill_gradientn(
    colors = c("red", "magenta", "blue"), # Your red-to-blue color scheme
    name = "p.adjust"
  ) +
  
  # 5. Customize labels and titles
  labs(
    title = "GO for unique SIRT1_KO genes",
    x = "Count",
    y = NULL # We don't need a y-axis title
  ) +
  
  # 6. Apply a clean theme and customize it to match the reference
  theme_minimal(base_size = 12) +
  theme(
    # Move the group labels (strips) to the outside and align them
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    
    # Remove the background from the group labels
    strip.background = element_blank(),
    
    # Adjust panel spacing
    panel.spacing = unit(0.5, "lines"),
    
    # General plot aesthetics
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines
    panel.grid.minor.x = element_blank()
  )

# --- Step 4: Display the Final Plot ---
print(p_final)


########## Get the number of genes that are found in unique pathways 
# 1. Get the names of the pathways that are unique to SIRT1_KO
sirt1_unique_pathway_names <- unique_pathways_df %>%
  dplyr::filter(Cluster == "SIRT1_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df <- as.data.frame(compare_go_results)

# 3. Filter the full results for SIRT1_KO AND its unique pathways
sirt1_unique_genes_df <- full_results_df %>%
  dplyr::filter(
    Cluster == "SIRT1_KO",
    Description %in% sirt1_unique_pathway_names
  )

# 4. Extract, split, and combine all the gene IDs from the 'geneID' column
# The geneID column is a single string with genes separated by "/"
all_sirt1_unique_genes <- sirt1_unique_genes_df %>%
  pull(geneID) %>%              # Extract the geneID column
  strsplit(split = "/") %>%     # Split each string by the "/"
  unlist() %>%                  # Flatten the list of lists into a single vector
  unique()                      # Keep only the unique gene IDs

# 5. Count them!
number_of_sirt1_genes_in_unique_pathways <- length(all_sirt1_unique_genes)

cat("The number of unique genes contributing to the unique functions of SIRT1 is:",
    number_of_sirt1_genes_in_unique_pathways, "\n")




######################## Extract the table for SIRT2
sirt2_unique_df <- list_of_unique_pathway_tables$SIRT2_KO

# Check if there's anything to process
if (!is.null(sirt2_unique_df) && nrow(sirt2_unique_df) > 0) {
  
  # Извлекаем список генов для SIRT1 из исходной таблицы результатов.
  # Это нужно для слота @gene в новом объекте.
  sirt2_gene_list <- gene_cluster_entrez$SIRT2_KO # Предполагается, что у вас есть этот список
  
  # Создаем новый объект класса enrichResult вручную.
  # Мы заполняем его данными из нашей таблицы sirt1_unique_df.
  sirt2_enrich_object_clean <- new("enrichResult",
                                   result = sirt2_unique_df,
                                   pvalueCutoff = 0.05, # Эти значения могут быть любыми
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.1,
                                   gene = as.character(sirt2_gene_list),
                                   universe = as.character(background_entrez_ids),
                                   organism = "HOMO SAPIENS",
                                   keytype = "ENTREZID",
                                   ontology = "BP"
  )
  
  # --- Шаг 4: Упрощаем (simplify) этот новый, чистый объект ---
  sirt2_simplified_results <- simplify(
    sirt2_enrich_object_clean,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
  # --- Шаг 5: Создаем barplot ---
  # Теперь barplot получит объект того класса, который он ожидает.
  
  p_sirt2_barplot <- barplot(
    sirt2_simplified_results,
    showCategory = 15,
    colorBy = "p.adjust",
    title = "GO for unique SIRT2_KO genes"
  ) +
    scale_fill_gradient(low = "red", high = "blue")
  
  # --- Шаг 6: Отображаем график ---
  print(p_sirt2_barplot)
  
} else {
  cat("Не найдено уникальных путей для SIRT2_KO для построения графика.\n")
}

# Преобразуем результаты в data.frame
sirt2_results_df <- as.data.frame(sirt2_simplified_results)

# Проверяем, что в таблице есть строки
if (nrow(sirt2_results_df) > 0) {
  
  # 3. Сохраняем таблицу в CSV файл
  # Этот файл вы будете редактировать вручную
  write.csv(sirt2_results_df, "SIRT2_unique_down_GO_for_grouping.csv", row.names = FALSE)
  
  cat("Таблица 'SIRT2_unique_down_GO_for_grouping.csv' была сохранена в вашей рабочей директории.\n")
  cat("Теперь откройте ее в Excel или Google Sheets для следующего шага.\n")
  
} else {
  cat("Нет данных для анализа и группировки.\n")
}


# --- 1. Читаем вашу новую, сгруппированную таблицу ---
grouped_data <- read.csv2("~/SIRT17_hMSCs/DESeq/SIRT2_unique_down_GO_grouped.csv")

# --- Step 2: Prepare the Data for Plotting ---

# First, ensure the numeric columns are correctly formatted
plot_data <- grouped_data %>%
  mutate(
    # Make sure p.adjust and Count are numeric
    p.adjust = as.numeric(gsub(",", ".", p.adjust)),
    Count = as.numeric(Count),
    # Filter out any rows that failed to parse
    !is.na(p.adjust),
    !is.na(Count)
  )

# For clarity, let's select a limited number of top pathways per theme
plot_data <- plot_data %>%
  group_by(Theme) %>%
  # Arrange by p.adjust within each group and take the top 5
  arrange(p.adjust, .by_group = TRUE) %>%
  slice(1:5) %>%
  ungroup()

# --- CRITICAL STEP for Ordering ---
# To plot correctly, we need to order the 'Description' factor levels
# based on their group (Theme) and their value (Count or p.adjust).
# We also need to order the 'Theme' factor levels.
# Let's order Themes by the max count within each theme for a nice visual.
theme_order <- plot_data %>%
  group_by(Theme) %>%
  summarise(max_count = max(Count)) %>%
  arrange(desc(max_count)) %>%
  pull(Theme)

plot_data <- plot_data %>%
  mutate(
    # Order the main pathway descriptions
    Description = reorder_within(Description, Count, Theme),
    # Order the Theme factor itself
    Theme = factor(Theme, levels = rev(theme_order))
  )


# --- Step 3: Create the Customized Bar Plot ---

p_final <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
  
  # 1. Create the bar plot
  geom_bar(stat = "identity") +
  
  # 2. Add the faceting by Theme. This creates the groups.
  # The 'switch = "y"' argument is the key to moving labels to the right.
  facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  # 3. Apply the correct sorting for the y-axis
  scale_y_reordered() +
  
  # 4. Customize the color scale
  scale_fill_gradientn(
    colors = c("red", "magenta", "blue"), # Your red-to-blue color scheme
    name = "p.adjust"
  ) +
  
  # 5. Customize labels and titles
  labs(
    title = "GO for unique SIRT2_KO genes",
    x = "Count",
    y = NULL # We don't need a y-axis title
  ) +
  
  # 6. Apply a clean theme and customize it to match the reference
  theme_minimal(base_size = 12) +
  theme(
    # Move the group labels (strips) to the outside and align them
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    
    # Remove the background from the group labels
    strip.background = element_blank(),
    
    # Adjust panel spacing
    panel.spacing = unit(0.5, "lines"),
    
    # General plot aesthetics
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines
    panel.grid.minor.x = element_blank()
  )

# --- Step 4: Display the Final Plot ---
print(p_final)


########## Get the number of genes that are found in unique pathways 
# 1. Get the names of the pathways that are unique to SIRT2_KO
sirt2_unique_pathway_names <- unique_pathways_df %>%
  dplyr::filter(Cluster == "SIRT2_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df <- as.data.frame(compare_go_results)

# 3. Filter the full results for SIRT2_KO AND its unique pathways
sirt2_unique_genes_df <- full_results_df %>%
  dplyr::filter(
    Cluster == "SIRT2_KO",
    Description %in% sirt2_unique_pathway_names
  )

# 4. Extract, split, and combine all the gene IDs from the 'geneID' column
# The geneID column is a single string with genes separated by "/"
all_sirt2_unique_genes <- sirt2_unique_genes_df %>%
  pull(geneID) %>%              # Extract the geneID column
  strsplit(split = "/") %>%     # Split each string by the "/"
  unlist() %>%                  # Flatten the list of lists into a single vector
  unique()                      # Keep only the unique gene IDs

# 5. Count them!
number_of_sirt2_genes_in_unique_pathways <- length(all_sirt2_unique_genes)

cat("The number of unique genes contributing to the unique functions of SIRT2 is:",
    number_of_sirt2_genes_in_unique_pathways, "\n")



######################## Extract the table for SIRT3
sirt3_unique_df <- list_of_unique_pathway_tables$SIRT3_KO

# Check if there's anything to process
if (!is.null(sirt3_unique_df) && nrow(sirt3_unique_df) > 0) {
  
  # Извлекаем список генов для SIRT1 из исходной таблицы результатов.
  # Это нужно для слота @gene в новом объекте.
  sirt3_gene_list <- gene_cluster_entrez$SIRT3_KO # Предполагается, что у вас есть этот список
  
  # Создаем новый объект класса enrichResult вручную.
  # Мы заполняем его данными из нашей таблицы sirt1_unique_df.
  sirt3_enrich_object_clean <- new("enrichResult",
                                   result = sirt3_unique_df,
                                   pvalueCutoff = 0.05, # Эти значения могут быть любыми
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.1,
                                   gene = as.character(sirt3_gene_list),
                                   universe = as.character(background_entrez_ids),
                                   organism = "HOMO SAPIENS",
                                   keytype = "ENTREZID",
                                   ontology = "BP"
  )
  
  # --- Шаг 4: Упрощаем (simplify) этот новый, чистый объект ---
  sirt3_simplified_results <- simplify(
    sirt3_enrich_object_clean,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
  # --- Шаг 5: Создаем barplot ---
  # Теперь barplot получит объект того класса, который он ожидает.
  
  p_sirt3_barplot <- barplot(
    sirt3_simplified_results,
    showCategory = 15,
    colorBy = "p.adjust",
    title = "GO for unique SIRT3_KO genes"
  ) +
    scale_fill_gradient(low = "red", high = "blue")
  
  # --- Шаг 6: Отображаем график ---
  print(p_sirt3_barplot)
  
} else {
  cat("Не найдено уникальных путей для SIRT3_KO для построения графика.\n")
}

# Преобразуем результаты в data.frame
sirt3_results_df <- as.data.frame(sirt3_simplified_results)

# Проверяем, что в таблице есть строки
if (nrow(sirt3_results_df) > 0) {
  
  # 3. Сохраняем таблицу в CSV файл
  # Этот файл вы будете редактировать вручную
  write.csv(sirt3_results_df, "SIRT3_unique_down_GO_for_grouping.csv", row.names = FALSE)
  
  cat("Таблица 'SIRT3_unique_down_GO_for_grouping.csv' была сохранена в вашей рабочей директории.\n")
  cat("Теперь откройте ее в Excel или Google Sheets для следующего шага.\n")
  
} else {
  cat("Нет данных для анализа и группировки.\n")
}


# --- 1. Читаем вашу новую, сгруппированную таблицу ---
grouped_data <- read.csv2("~/SIRT17_hMSCs/DESeq/SIRT3_unique_down_GO_grouped.csv")

# --- Step 2: Prepare the Data for Plotting ---

# First, ensure the numeric columns are correctly formatted
plot_data <- grouped_data %>%
  mutate(
    # Make sure p.adjust and Count are numeric
    p.adjust = as.numeric(gsub(",", ".", p.adjust)),
    Count = as.numeric(Count),
    # Filter out any rows that failed to parse
    !is.na(p.adjust),
    !is.na(Count)
  )

# For clarity, let's select a limited number of top pathways per theme
plot_data <- plot_data %>%
  group_by(Theme) %>%
  # Arrange by p.adjust within each group and take the top 5
  arrange(p.adjust, .by_group = TRUE) %>%
  slice(1:5) %>%
  ungroup()

# --- CRITICAL STEP for Ordering ---
# To plot correctly, we need to order the 'Description' factor levels
# based on their group (Theme) and their value (Count or p.adjust).
# We also need to order the 'Theme' factor levels.
# Let's order Themes by the max count within each theme for a nice visual.
theme_order <- plot_data %>%
  group_by(Theme) %>%
  summarise(max_count = max(Count)) %>%
  arrange(desc(max_count)) %>%
  pull(Theme)

plot_data <- plot_data %>%
  mutate(
    # Order the main pathway descriptions
    Description = reorder_within(Description, Count, Theme),
    # Order the Theme factor itself
    Theme = factor(Theme, levels = rev(theme_order))
  )


# --- Step 3: Create the Customized Bar Plot ---

p_final <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
  
  # 1. Create the bar plot
  geom_bar(stat = "identity") +
  
  # 2. Add the faceting by Theme. This creates the groups.
  # The 'switch = "y"' argument is the key to moving labels to the right.
  facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  # 3. Apply the correct sorting for the y-axis
  scale_y_reordered() +
  
  # 4. Customize the color scale
  scale_fill_gradientn(
    colors = c("red", "magenta", "blue"), # Your red-to-blue color scheme
    name = "p.adjust"
  ) +
  
  # 5. Customize labels and titles
  labs(
    title = "GO for unique SIRT3_KO genes",
    x = "Count",
    y = NULL # We don't need a y-axis title
  ) +
  
  # 6. Apply a clean theme and customize it to match the reference
  theme_minimal(base_size = 12) +
  theme(
    # Move the group labels (strips) to the outside and align them
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    
    # Remove the background from the group labels
    strip.background = element_blank(),
    
    # Adjust panel spacing
    panel.spacing = unit(0.5, "lines"),
    
    # General plot aesthetics
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines
    panel.grid.minor.x = element_blank()
  )

# --- Step 4: Display the Final Plot ---
print(p_final)



########## Get the number of genes that are found in unique pathways 
# 1. Get the names of the pathways that are unique to SIRT3_KO
sirt3_unique_pathway_names <- unique_pathways_df %>%
  dplyr::filter(Cluster == "SIRT3_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df <- as.data.frame(compare_go_results)

# 3. Filter the full results for SIRT3_KO AND its unique pathways
sirt3_unique_genes_df <- full_results_df %>%
  dplyr::filter(
    Cluster == "SIRT3_KO",
    Description %in% sirt3_unique_pathway_names
  )

# 4. Extract, split, and combine all the gene IDs from the 'geneID' column
# The geneID column is a single string with genes separated by "/"
all_sirt3_unique_genes <- sirt3_unique_genes_df %>%
  pull(geneID) %>%              # Extract the geneID column
  strsplit(split = "/") %>%     # Split each string by the "/"
  unlist() %>%                  # Flatten the list of lists into a single vector
  unique()                      # Keep only the unique gene IDs

# 5. Count them!
number_of_sirt3_genes_in_unique_pathways <- length(all_sirt3_unique_genes)

cat("The number of unique genes contributing to the unique functions of SIRT3 is:",
    number_of_sirt3_genes_in_unique_pathways, "\n")



######################## Extract the table for SIRT4
sirt4_unique_df <- list_of_unique_pathway_tables$SIRT4_KO

# Check if there's anything to process
if (!is.null(sirt4_unique_df) && nrow(sirt4_unique_df) > 0) {
  
  # Извлекаем список генов для SIRT1 из исходной таблицы результатов.
  # Это нужно для слота @gene в новом объекте.
  sirt4_gene_list <- gene_cluster_entrez$SIRT4_KO # Предполагается, что у вас есть этот список
  
  # Создаем новый объект класса enrichResult вручную.
  # Мы заполняем его данными из нашей таблицы sirt1_unique_df.
  sirt4_enrich_object_clean <- new("enrichResult",
                                   result = sirt4_unique_df,
                                   pvalueCutoff = 0.05, # Эти значения могут быть любыми
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.1,
                                   gene = as.character(sirt4_gene_list),
                                   universe = as.character(background_entrez_ids),
                                   organism = "HOMO SAPIENS",
                                   keytype = "ENTREZID",
                                   ontology = "BP"
  )
  
  # --- Шаг 4: Упрощаем (simplify) этот новый, чистый объект ---
  sirt4_simplified_results <- simplify(
    sirt4_enrich_object_clean,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
  # --- Шаг 5: Создаем barplot ---
  # Теперь barplot получит объект того класса, который он ожидает.
  
  p_sirt4_barplot <- barplot(
    sirt4_simplified_results,
    showCategory = 15,
    colorBy = "p.adjust",
    title = "GO for unique SIRT4_KO genes"
  ) +
    scale_fill_gradient(low = "red", high = "blue")
  
  # --- Шаг 6: Отображаем график ---
  print(p_sirt4_barplot)
  
} else {
  cat("Не найдено уникальных путей для SIRT4_KO для построения графика.\n")
}

# Преобразуем результаты в data.frame
sirt4_results_df <- as.data.frame(sirt4_simplified_results)

# Проверяем, что в таблице есть строки
if (nrow(sirt4_results_df) > 0) {
  
  # 3. Сохраняем таблицу в CSV файл
  # Этот файл вы будете редактировать вручную
  write.csv(sirt4_results_df, "SIRT4_unique_down_GO_for_grouping.csv", row.names = FALSE)
  
  cat("Таблица 'SIRT4_unique_down_GO_for_grouping.csv' была сохранена в вашей рабочей директории.\n")
  cat("Теперь откройте ее в Excel или Google Sheets для следующего шага.\n")
  
} else {
  cat("Нет данных для анализа и группировки.\n")
}


# --- 1. Читаем вашу новую, сгруппированную таблицу ---
grouped_data <- read.csv2("~/SIRT17_hMSCs/DESeq/SIRT4_unique_down_GO_grouped.csv")

# --- Step 2: Prepare the Data for Plotting ---

# First, ensure the numeric columns are correctly formatted
plot_data <- grouped_data %>%
  mutate(
    # Make sure p.adjust and Count are numeric
    p.adjust = as.numeric(gsub(",", ".", p.adjust)),
    Count = as.numeric(Count),
    # Filter out any rows that failed to parse
    !is.na(p.adjust),
    !is.na(Count)
  )

# For clarity, let's select a limited number of top pathways per theme
plot_data <- plot_data %>%
  group_by(Theme) %>%
  # Arrange by p.adjust within each group and take the top 5
  arrange(p.adjust, .by_group = TRUE) %>%
  slice(1:5) %>%
  ungroup()

# --- CRITICAL STEP for Ordering ---
# To plot correctly, we need to order the 'Description' factor levels
# based on their group (Theme) and their value (Count or p.adjust).
# We also need to order the 'Theme' factor levels.
# Let's order Themes by the max count within each theme for a nice visual.
theme_order <- plot_data %>%
  group_by(Theme) %>%
  summarise(max_count = max(Count)) %>%
  arrange(desc(max_count)) %>%
  pull(Theme)

plot_data <- plot_data %>%
  mutate(
    # Order the main pathway descriptions
    Description = reorder_within(Description, Count, Theme),
    # Order the Theme factor itself
    Theme = factor(Theme, levels = rev(theme_order))
  )


# --- Step 3: Create the Customized Bar Plot ---

p_final <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
  
  # 1. Create the bar plot
  geom_bar(stat = "identity") +
  
  # 2. Add the faceting by Theme. This creates the groups.
  # The 'switch = "y"' argument is the key to moving labels to the right.
  facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  # 3. Apply the correct sorting for the y-axis
  scale_y_reordered() +
  
  # 4. Customize the color scale
  scale_fill_gradientn(
    colors = c("red", "magenta", "blue"), # Your red-to-blue color scheme
    name = "p.adjust"
  ) +
  
  # 5. Customize labels and titles
  labs(
    title = "GO for unique SIRT4_KO genes",
    x = "Count",
    y = NULL # We don't need a y-axis title
  ) +
  
  # 6. Apply a clean theme and customize it to match the reference
  theme_minimal(base_size = 12) +
  theme(
    # Move the group labels (strips) to the outside and align them
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    
    # Remove the background from the group labels
    strip.background = element_blank(),
    
    # Adjust panel spacing
    panel.spacing = unit(0.5, "lines"),
    
    # General plot aesthetics
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines
    panel.grid.minor.x = element_blank()
  )

# --- Step 4: Display the Final Plot ---
print(p_final)



########## Get the number of genes that are found in unique pathways 
# 1. Get the names of the pathways that are unique to SIRT4_KO
sirt4_unique_pathway_names <- unique_pathways_df %>%
  dplyr::filter(Cluster == "SIRT4_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df <- as.data.frame(compare_go_results)

# 3. Filter the full results for SIRT4_KO AND its unique pathways
sirt4_unique_genes_df <- full_results_df %>%
  dplyr::filter(
    Cluster == "SIRT4_KO",
    Description %in% sirt4_unique_pathway_names
  )

# 4. Extract, split, and combine all the gene IDs from the 'geneID' column
# The geneID column is a single string with genes separated by "/"
all_sirt4_unique_genes <- sirt4_unique_genes_df %>%
  pull(geneID) %>%              # Extract the geneID column
  strsplit(split = "/") %>%     # Split each string by the "/"
  unlist() %>%                  # Flatten the list of lists into a single vector
  unique()                      # Keep only the unique gene IDs

# 5. Count them!
number_of_sirt4_genes_in_unique_pathways <- length(all_sirt4_unique_genes)

cat("The number of unique genes contributing to the unique functions of SIRT4 is:",
    number_of_sirt4_genes_in_unique_pathways, "\n")



######################## Extract the table for SIRT5
sirt5_unique_df <- list_of_unique_pathway_tables$SIRT5_KO

# Check if there's anything to process
if (!is.null(sirt5_unique_df) && nrow(sirt5_unique_df) > 0) {
  
  # Извлекаем список генов для SIRT1 из исходной таблицы результатов.
  # Это нужно для слота @gene в новом объекте.
  sirt5_gene_list <- gene_cluster_entrez$SIRT5_KO # Предполагается, что у вас есть этот список
  
  # Создаем новый объект класса enrichResult вручную.
  # Мы заполняем его данными из нашей таблицы sirt1_unique_df.
  sirt5_enrich_object_clean <- new("enrichResult",
                                   result = sirt5_unique_df,
                                   pvalueCutoff = 0.05, # Эти значения могут быть любыми
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.1,
                                   gene = as.character(sirt5_gene_list),
                                   universe = as.character(background_entrez_ids),
                                   organism = "HOMO SAPIENS",
                                   keytype = "ENTREZID",
                                   ontology = "BP"
  )
  
  # --- Шаг 4: Упрощаем (simplify) этот новый, чистый объект ---
  sirt5_simplified_results <- simplify(
    sirt5_enrich_object_clean,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
  # --- Шаг 5: Создаем barplot ---
  # Теперь barplot получит объект того класса, который он ожидает.
  
  p_sirt5_barplot <- barplot(
    sirt5_simplified_results,
    showCategory = 15,
    colorBy = "p.adjust",
    title = "GO for unique SIRT5_KO genes"
  ) +
    scale_fill_gradient(low = "red", high = "blue")
  
  # --- Шаг 6: Отображаем график ---
  print(p_sirt5_barplot)
  
} else {
  cat("Не найдено уникальных путей для SIRT5_KO для построения графика.\n")
}

# Преобразуем результаты в data.frame
sirt5_results_df <- as.data.frame(sirt5_simplified_results)

# Проверяем, что в таблице есть строки
if (nrow(sirt5_results_df) > 0) {
  
  # 3. Сохраняем таблицу в CSV файл
  # Этот файл вы будете редактировать вручную
  write.csv(sirt5_results_df, "SIRT5_unique_down_GO_for_grouping.csv", row.names = FALSE)
  
  cat("Таблица 'SIRT5_unique_down_GO_for_grouping.csv' была сохранена в вашей рабочей директории.\n")
  cat("Теперь откройте ее в Excel или Google Sheets для следующего шага.\n")
  
} else {
  cat("Нет данных для анализа и группировки.\n")
}


# --- 1. Читаем вашу новую, сгруппированную таблицу ---
grouped_data <- read.csv2("~/SIRT17_hMSCs/DESeq/SIRT5_unique_down_GO_grouped.csv")

# --- Step 2: Prepare the Data for Plotting ---

# First, ensure the numeric columns are correctly formatted
plot_data <- grouped_data %>%
  mutate(
    # Make sure p.adjust and Count are numeric
    p.adjust = as.numeric(gsub(",", ".", p.adjust)),
    Count = as.numeric(Count),
    # Filter out any rows that failed to parse
    !is.na(p.adjust),
    !is.na(Count)
  )

# For clarity, let's select a limited number of top pathways per theme
plot_data <- plot_data %>%
  group_by(Theme) %>%
  # Arrange by p.adjust within each group and take the top 5
  arrange(p.adjust, .by_group = TRUE) %>%
  slice(1:5) %>%
  ungroup()

# --- CRITICAL STEP for Ordering ---
# To plot correctly, we need to order the 'Description' factor levels
# based on their group (Theme) and their value (Count or p.adjust).
# We also need to order the 'Theme' factor levels.
# Let's order Themes by the max count within each theme for a nice visual.
theme_order <- plot_data %>%
  group_by(Theme) %>%
  summarise(max_count = max(Count)) %>%
  arrange(desc(max_count)) %>%
  pull(Theme)

plot_data <- plot_data %>%
  mutate(
    # Order the main pathway descriptions
    Description = reorder_within(Description, Count, Theme),
    # Order the Theme factor itself
    Theme = factor(Theme, levels = rev(theme_order))
  )


# --- Step 3: Create the Customized Bar Plot ---

p_final <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
  
  # 1. Create the bar plot
  geom_bar(stat = "identity") +
  
  # 2. Add the faceting by Theme. This creates the groups.
  # The 'switch = "y"' argument is the key to moving labels to the right.
  facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  # 3. Apply the correct sorting for the y-axis
  scale_y_reordered() +
  
  # 4. Customize the color scale
  scale_fill_gradientn(
    colors = c("red", "magenta", "blue"), # Your red-to-blue color scheme
    name = "p.adjust"
  ) +
  
  # 5. Customize labels and titles
  labs(
    title = "GO for unique SIRT5_KO genes",
    x = "Count",
    y = NULL # We don't need a y-axis title
  ) +
  
  # 6. Apply a clean theme and customize it to match the reference
  theme_minimal(base_size = 12) +
  theme(
    # Move the group labels (strips) to the outside and align them
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    
    # Remove the background from the group labels
    strip.background = element_blank(),
    
    # Adjust panel spacing
    panel.spacing = unit(0.5, "lines"),
    
    # General plot aesthetics
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines
    panel.grid.minor.x = element_blank()
  )

# --- Step 4: Display the Final Plot ---
print(p_final)



########## Get the number of genes that are found in unique pathways 
# 1. Get the names of the pathways that are unique to SIRT5_KO
sirt5_unique_pathway_names <- unique_pathways_df %>%
  dplyr::filter(Cluster == "SIRT5_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df <- as.data.frame(compare_go_results)

# 3. Filter the full results for SIRT5_KO AND its unique pathways
sirt5_unique_genes_df <- full_results_df %>%
  dplyr::filter(
    Cluster == "SIRT5_KO",
    Description %in% sirt5_unique_pathway_names
  )

# 4. Extract, split, and combine all the gene IDs from the 'geneID' column
# The geneID column is a single string with genes separated by "/"
all_sirt5_unique_genes <- sirt5_unique_genes_df %>%
  pull(geneID) %>%              # Extract the geneID column
  strsplit(split = "/") %>%     # Split each string by the "/"
  unlist() %>%                  # Flatten the list of lists into a single vector
  unique()                      # Keep only the unique gene IDs

# 5. Count them!
number_of_sirt5_genes_in_unique_pathways <- length(all_sirt5_unique_genes)

cat("The number of unique genes contributing to the unique functions of SIRT5 is:",
    number_of_sirt5_genes_in_unique_pathways, "\n")



######################## Extract the table for SIRT6
sirt6_unique_df <- list_of_unique_pathway_tables$SIRT6_KO

# Check if there's anything to process
if (!is.null(sirt6_unique_df) && nrow(sirt6_unique_df) > 0) {
  
  # Извлекаем список генов для SIRT1 из исходной таблицы результатов.
  # Это нужно для слота @gene в новом объекте.
  sirt6_gene_list <- gene_cluster_entrez$SIRT6_KO # Предполагается, что у вас есть этот список
  
  # Создаем новый объект класса enrichResult вручную.
  # Мы заполняем его данными из нашей таблицы sirt1_unique_df.
  sirt6_enrich_object_clean <- new("enrichResult",
                                   result = sirt6_unique_df,
                                   pvalueCutoff = 0.05, # Эти значения могут быть любыми
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.1,
                                   gene = as.character(sirt6_gene_list),
                                   universe = as.character(background_entrez_ids),
                                   organism = "HOMO SAPIENS",
                                   keytype = "ENTREZID",
                                   ontology = "BP"
  )
  
  # --- Шаг 4: Упрощаем (simplify) этот новый, чистый объект ---
  sirt6_simplified_results <- simplify(
    sirt6_enrich_object_clean,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
  # --- Шаг 5: Создаем barplot ---
  # Теперь barplot получит объект того класса, который он ожидает.
  
  p_sirt6_barplot <- barplot(
    sirt6_simplified_results,
    showCategory = 15,
    colorBy = "p.adjust",
    title = "GO for unique SIRT6_KO genes"
  ) +
    scale_fill_gradient(low = "red", high = "blue")
  
  # --- Шаг 6: Отображаем график ---
  print(p_sirt6_barplot)
  
} else {
  cat("Не найдено уникальных путей для SIRT6_KO для построения графика.\n")
}

# Преобразуем результаты в data.frame
sirt6_results_df <- as.data.frame(sirt6_simplified_results)

# Проверяем, что в таблице есть строки
if (nrow(sirt6_results_df) > 0) {
  
  # 3. Сохраняем таблицу в CSV файл
  # Этот файл вы будете редактировать вручную
  write.csv(sirt6_results_df, "SIRT6_unique_down_GO_for_grouping.csv", row.names = FALSE)
  
  cat("Таблица 'SIRT6_unique_down_GO_for_grouping.csv' была сохранена в вашей рабочей директории.\n")
  cat("Теперь откройте ее в Excel или Google Sheets для следующего шага.\n")
  
} else {
  cat("Нет данных для анализа и группировки.\n")
}


# --- 1. Читаем вашу новую, сгруппированную таблицу ---
grouped_data <- read.csv2("~/SIRT17_hMSCs/DESeq/SIRT6_unique_down_GO_grouped.csv")

# --- Step 2: Prepare the Data for Plotting ---

# First, ensure the numeric columns are correctly formatted
plot_data <- grouped_data %>%
  mutate(
    # Make sure p.adjust and Count are numeric
    p.adjust = as.numeric(gsub(",", ".", p.adjust)),
    Count = as.numeric(Count),
    # Filter out any rows that failed to parse
    !is.na(p.adjust),
    !is.na(Count)
  )

# For clarity, let's select a limited number of top pathways per theme
plot_data <- plot_data %>%
  group_by(Theme) %>%
  # Arrange by p.adjust within each group and take the top 5
  arrange(p.adjust, .by_group = TRUE) %>%
  slice(1:5) %>%
  ungroup()

# --- CRITICAL STEP for Ordering ---
# To plot correctly, we need to order the 'Description' factor levels
# based on their group (Theme) and their value (Count or p.adjust).
# We also need to order the 'Theme' factor levels.
# Let's order Themes by the max count within each theme for a nice visual.
theme_order <- plot_data %>%
  group_by(Theme) %>%
  summarise(max_count = max(Count)) %>%
  arrange(desc(max_count)) %>%
  pull(Theme)

plot_data <- plot_data %>%
  mutate(
    # Order the main pathway descriptions
    Description = reorder_within(Description, Count, Theme),
    # Order the Theme factor itself
    Theme = factor(Theme, levels = rev(theme_order))
  )


# --- Step 3: Create the Customized Bar Plot ---

p_final <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
  
  # 1. Create the bar plot
  geom_bar(stat = "identity") +
  
  # 2. Add the faceting by Theme. This creates the groups.
  # The 'switch = "y"' argument is the key to moving labels to the right.
  facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  # 3. Apply the correct sorting for the y-axis
  scale_y_reordered() +
  
  # 4. Customize the color scale
  scale_fill_gradientn(
    colors = c("red", "magenta", "blue"), # Your red-to-blue color scheme
    name = "p.adjust"
  ) +
  
  # 5. Customize labels and titles
  labs(
    title = "GO for unique SIRT6_KO genes",
    x = "Count",
    y = NULL # We don't need a y-axis title
  ) +
  
  # 6. Apply a clean theme and customize it to match the reference
  theme_minimal(base_size = 12) +
  theme(
    # Move the group labels (strips) to the outside and align them
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    
    # Remove the background from the group labels
    strip.background = element_blank(),
    
    # Adjust panel spacing
    panel.spacing = unit(0.5, "lines"),
    
    # General plot aesthetics
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines
    panel.grid.minor.x = element_blank()
  )

# --- Step 4: Display the Final Plot ---
print(p_final)



########## Get the number of genes that are found in unique pathways 
# 1. Get the names of the pathways that are unique to SIRT6_KO
sirt6_unique_pathway_names <- unique_pathways_df %>%
  dplyr::filter(Cluster == "SIRT6_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df <- as.data.frame(compare_go_results)

# 3. Filter the full results for SIRT6_KO AND its unique pathways
sirt6_unique_genes_df <- full_results_df %>%
  dplyr::filter(
    Cluster == "SIRT6_KO",
    Description %in% sirt6_unique_pathway_names
  )

# 4. Extract, split, and combine all the gene IDs from the 'geneID' column
# The geneID column is a single string with genes separated by "/"
all_sirt6_unique_genes <- sirt6_unique_genes_df %>%
  pull(geneID) %>%              # Extract the geneID column
  strsplit(split = "/") %>%     # Split each string by the "/"
  unlist() %>%                  # Flatten the list of lists into a single vector
  unique()                      # Keep only the unique gene IDs

# 5. Count them!
number_of_sirt6_genes_in_unique_pathways <- length(all_sirt6_unique_genes)

cat("The number of unique genes contributing to the unique functions of SIRT6 is:",
    number_of_sirt6_genes_in_unique_pathways, "\n")




######################## Extract the table for SIRT7
sirt7_unique_df <- list_of_unique_pathway_tables$SIRT7_KO

# Check if there's anything to process
if (!is.null(sirt7_unique_df) && nrow(sirt7_unique_df) > 0) {
  
  # Извлекаем список генов для SIRT1 из исходной таблицы результатов.
  # Это нужно для слота @gene в новом объекте.
  sirt7_gene_list <- gene_cluster_entrez$SIRT7_KO # Предполагается, что у вас есть этот список
  
  # Создаем новый объект класса enrichResult вручную.
  # Мы заполняем его данными из нашей таблицы sirt1_unique_df.
  sirt7_enrich_object_clean <- new("enrichResult",
                                   result = sirt7_unique_df,
                                   pvalueCutoff = 0.05, # Эти значения могут быть любыми
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.1,
                                   gene = as.character(sirt7_gene_list),
                                   universe = as.character(background_entrez_ids),
                                   organism = "HOMO SAPIENS",
                                   keytype = "ENTREZID",
                                   ontology = "BP"
  )
  
  # --- Шаг 4: Упрощаем (simplify) этот новый, чистый объект ---
  sirt7_simplified_results <- simplify(
    sirt7_enrich_object_clean,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
  # --- Шаг 5: Создаем barplot ---
  # Теперь barplot получит объект того класса, который он ожидает.
  
  p_sirt7_barplot <- barplot(
    sirt7_simplified_results,
    showCategory = 15,
    colorBy = "p.adjust",
    title = "GO for unique SIRT7_KO genes"
  ) +
    scale_fill_gradient(low = "red", high = "blue")
  
  # --- Шаг 6: Отображаем график ---
  print(p_sirt7_barplot)
  
} else {
  cat("Не найдено уникальных путей для SIRT7_KO для построения графика.\n")
}

# Преобразуем результаты в data.frame
sirt7_results_df <- as.data.frame(sirt7_simplified_results)

# Проверяем, что в таблице есть строки
if (nrow(sirt7_results_df) > 0) {
  
  # 3. Сохраняем таблицу в CSV файл
  # Этот файл вы будете редактировать вручную
  write.csv(sirt7_results_df, "SIRT7_unique_down_GO_for_grouping.csv", row.names = FALSE)
  
  cat("Таблица 'SIRT7_unique_down_GO_for_grouping.csv' была сохранена в вашей рабочей директории.\n")
  cat("Теперь откройте ее в Excel или Google Sheets для следующего шага.\n")
  
} else {
  cat("Нет данных для анализа и группировки.\n")
}


# --- 1. Читаем вашу новую, сгруппированную таблицу ---
grouped_data <- read.csv2("~/SIRT17_hMSCs/DESeq/SIRT7_unique_down_GO_grouped.csv")

# --- Step 2: Prepare the Data for Plotting ---

# First, ensure the numeric columns are correctly formatted
plot_data <- grouped_data %>%
  mutate(
    # Make sure p.adjust and Count are numeric
    p.adjust = as.numeric(gsub(",", ".", p.adjust)),
    Count = as.numeric(Count),
    # Filter out any rows that failed to parse
    !is.na(p.adjust),
    !is.na(Count)
  )

# For clarity, let's select a limited number of top pathways per theme
plot_data <- plot_data %>%
  group_by(Theme) %>%
  # Arrange by p.adjust within each group and take the top 5
  arrange(p.adjust, .by_group = TRUE) %>%
  slice(1:5) %>%
  ungroup()

# --- CRITICAL STEP for Ordering ---
# To plot correctly, we need to order the 'Description' factor levels
# based on their group (Theme) and their value (Count or p.adjust).
# We also need to order the 'Theme' factor levels.
# Let's order Themes by the max count within each theme for a nice visual.
theme_order <- plot_data %>%
  group_by(Theme) %>%
  summarise(max_count = max(Count)) %>%
  arrange(desc(max_count)) %>%
  pull(Theme)

plot_data <- plot_data %>%
  mutate(
    # Order the main pathway descriptions
    Description = reorder_within(Description, Count, Theme),
    # Order the Theme factor itself
    Theme = factor(Theme, levels = rev(theme_order))
  )


# --- Step 3: Create the Customized Bar Plot ---

p_final <- ggplot(plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
  
  # 1. Create the bar plot
  geom_bar(stat = "identity") +
  
  # 2. Add the faceting by Theme. This creates the groups.
  # The 'switch = "y"' argument is the key to moving labels to the right.
  facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  # 3. Apply the correct sorting for the y-axis
  scale_y_reordered() +
  
  # 4. Customize the color scale
  scale_fill_gradientn(
    colors = c("red", "magenta", "blue"), # Your red-to-blue color scheme
    name = "p.adjust"
  ) +
  
  # 5. Customize labels and titles
  labs(
    title = "GO for unique SIRT6_KO genes",
    x = "Count",
    y = NULL # We don't need a y-axis title
  ) +
  
  # 6. Apply a clean theme and customize it to match the reference
  theme_minimal(base_size = 12) +
  theme(
    # Move the group labels (strips) to the outside and align them
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
    
    # Remove the background from the group labels
    strip.background = element_blank(),
    
    # Adjust panel spacing
    panel.spacing = unit(0.5, "lines"),
    
    # General plot aesthetics
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(), # Remove horizontal grid lines
    panel.grid.minor.x = element_blank()
  )

# --- Step 4: Display the Final Plot ---
print(p_final)



########## Get the number of genes that are found in unique pathways 
# 1. Get the names of the pathways that are unique to SIRT7_KO
sirt7_unique_pathway_names <- unique_pathways_df %>%
  dplyr::filter(Cluster == "SIRT7_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df <- as.data.frame(compare_go_results)

# 3. Filter the full results for SIRT7_KO AND its unique pathways
sirt7_unique_genes_df <- full_results_df %>%
  dplyr::filter(
    Cluster == "SIRT7_KO",
    Description %in% sirt7_unique_pathway_names
  )

# 4. Extract, split, and combine all the gene IDs from the 'geneID' column
# The geneID column is a single string with genes separated by "/"
all_sirt7_unique_genes <- sirt7_unique_genes_df %>%
  pull(geneID) %>%              # Extract the geneID column
  strsplit(split = "/") %>%     # Split each string by the "/"
  unlist() %>%                  # Flatten the list of lists into a single vector
  unique()                      # Keep only the unique gene IDs

# 5. Count them!
number_of_sirt7_genes_in_unique_pathways <- length(all_sirt7_unique_genes)

cat("The number of unique genes contributing to the unique functions of SIRT7 is:",
    number_of_sirt7_genes_in_unique_pathways, "\n")






#### VI. Analyze unique down-regul. SIRT3-5 (mitochondrial SIRT) genes against DNA repair pathways from MSigDB ####

# Load required libraries
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)

#### 1. SIRT3 #####
# unique_genes_full_list$SIRT3_KO - unique SIRT3 genes 
# downregulated_sig_lists$SIRT3 - down-regulated sig. SIRT3 DEGs

# Filter for down-regulated genes (negative log2FoldChange)
sirt3_down_genes <- unique_genes_full_list$SIRT3_KO %>% 
  filter(log2FoldChange < 0) %>% 
  pull(gene_name)

sirt3_down_g <- downregulated_sig_lists$SIRT3 %>%
  pull(gene_name)

# Convert SIRT3 down-regulated genes to Entrez IDs
sirt3_down_entrez_ids <- bitr(sirt3_down_genes,
                          fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)$ENTREZID

length(sirt3_down_entrez_ids) # 23

sirt3_down_g_entrez_ids <- bitr(
  sirt3_down_g,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)$ENTREZID

# Target Entrez IDs to extract (These are genes from DNA damage response unique GO pathways)
target_entrez <- c("641", "675", "2625", "4436", "64710", "7056", "7516", "7518")

# Extract only these genes from your vector
sirt3_extracted_entrez <- sirt3_down_g_entrez_ids[sirt3_down_g_entrez_ids %in% target_entrez]

# --- 2. Get DNA repair gene sets from MSigDB ---
# Get all human gene sets from Hallmark and C2 (Reactome) collections
msig_h <- msigdbr(species = "Homo sapiens", collection = "H")
msig_c2 <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")

# Check the column names to see what's available
print("MSigDB column names:")
print(colnames(msig_h)) # we need ncbi_gene column that contains ENTREZ ID 

# Filter for DNA repair-related pathways
dna_repair_terms <- c(
  # Hallmark DNA Repair
  "HALLMARK_DNA_REPAIR",
  
  # Reactome DNA Repair Pathways
  "REACTOME_HOMOLOGY_DIRECTED_REPAIR",
  "REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR",
  "REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ",
  "REACTOME_BASE_EXCISION_REPAIR",
  "REACTOME_NUCLEOTIDE_EXCISION_REPAIR",
  "REACTOME_MISMATCH_REPAIR",
  "REACTOME_GLOBAL_GENOME_NUCLEOTIDE_EXCISION_REPAIR_GG_NER",
  "REACTOME_TRANSCRIPTION_COUPLED_NUCLEOTIDE_EXCISION_REPAIR_TC_NER",
  "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR",
  "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS"
)

# Create the TERM2GENE data frame directly
term2gene <- data.frame()

for(term in dna_repair_terms) {
  if(term %in% unique(msig_h$gs_name)) {
    genes_df <- msig_h %>% 
      filter(gs_name == term) %>%
      dplyr::select(term = gs_name, gene = ncbi_gene)
  } else if(term %in% unique(msig_c2$gs_name)) {
    genes_df <- msig_c2 %>% 
      filter(gs_name == term) %>%
      dplyr::select(term = gs_name, gene = ncbi_gene)
  } else {
    warning(paste("Term not found:", term))
    next
  }
  term2gene <- rbind(term2gene, genes_df)
}

# Convert gene column to character (important for matching)
term2gene$gene <- as.character(term2gene$gene)

print("DNA repair gene sets summary:")
print(table(term2gene$term))

# --- 3. Check gene ID compatibility ---
# Check if any of our genes are in the DNA repair sets
intersection <- intersect(sirt3_extracted_entrez, term2gene$gene)
print(paste("Number of SIRT3 down genes in DNA repair pathways:", length(intersection)))

if(length(intersection) == 0) {
  print("No overlap found. Checking gene ID formats...")
  print(paste("Sample of SIRT3 Entrez IDs:", paste(head(sirt3_down_entrez_ids), collapse = ", ")))
  print(paste("Sample of DNA repair Entrez IDs:", paste(head(unique(term2gene$gene)), collapse = ", ")))
}

# --- 4. Perform enrichment analysis ---
if(length(intersection) > 0) {
  sirt3_dna_repair_enrich <- enricher(
    gene = sirt3_extracted_entrez,
    universe = background_entrez_ids,
    TERM2GENE = term2gene,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 3,    # Reduced for more sensitivity
    maxGSSize = 500
  )
  
  # --- 5. View and interpret results ---
  if(!is.null(sirt3_dna_repair_enrich)) {
    
    results_df <- sirt3_dna_repair_enrich@result
    significant_results <- results_df %>% filter(p.adjust < 1)
    
    print("SIRT3 Down-regulated Genes - DNA Repair Pathway Enrichment:")
    print(paste("Total pathways tested:", nrow(results_df)))
    print(paste("Significant pathways (p.adj < 0.05):", nrow(significant_results)))
    
    if(nrow(significant_results) > 0) {
      # Print significant results
      print("Significant DNA Repair Pathways:")
      print(significant_results %>% 
              dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count) %>%
              arrange(p.adjust))
      
      # Create visualization
      plot_data <- significant_results %>%
        arrange(p.adjust) %>%
        head(15)
      
      # Create dot plot - color by p.adjust directly
      p <- ggplot(plot_data, aes(x = Count, y = reorder(Description, Count), 
                                 size = Count, color = p.adjust)) +
        geom_point(alpha = 0.8) +
        scale_color_gradient(low = "red", high = "blue", 
                             name = "Adjusted p-value",
                             trans = "reverse") +  # Reverse so red = more significant
        scale_size_continuous(range = c(3, 8), name = "Gene Count") +
        labs(title = "SIRT3: DNA repair pathways enriched in down-regulated genes",
             subtitle = "Pathways normally activated by SIRT3",
             x = "Number of Genes", y = "") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5),
              axis.text.y = element_text(size = 10))
      
      print(p)
      
      # Extract and examine specific genes
      print("Genes driving significant enrichments:")
      for(i in 1:nrow(plot_data)) {
        pathway <- plot_data$ID[i]
        genes_in_pathway <- strsplit(plot_data$geneID[i], "/")[[1]]
        
        # Convert back to symbols for readability
        gene_symbols <- bitr(genes_in_pathway,
                             fromType = "ENTREZID", 
                             toType = "SYMBOL",
                             OrgDb = org.Hs.eg.db)$SYMBOL
        
        print(paste("=== ", pathway, " ==="))
        print(paste("Genes (", length(gene_symbols), "):", paste(gene_symbols, collapse = ", ")))
        cat("\n")
      }
      
    } else {
      print("No significant DNA repair pathways found.")
      print("Top pathways (even if not significant):")
      print(results_df %>% 
              arrange(p.adjust) %>% 
              head(10) %>%
              dplyr::select(ID, Description, p.adjust, Count))
    }
  }
} else {
  print("No overlap between SIRT3 down-regulated genes and DNA repair pathways.")
}



# Volcano plot for SIRT3
library(EnhancedVolcano)

# Filter to get only the 5 DNA repair genes
dna_repair_genes <- downregulated_sig_lists$SIRT3 %>%
  filter(gene_name %in% c('BLM', 'BRCA2', 'XRCC2', 'XRCC4', 'MSH2'))

# Simple volcano plot with just these genes
p3 <- EnhancedVolcano(dna_repair_genes,
                      lab = dna_repair_genes$gene_name, 
                      x = 'log2FoldChange',
                      y = 'padj', 
                      title = 'Down-regulated DNA repair genes in SIRT3 KO',
                      subtitle = "FC cutoff = 0.58; p-value cutoff = 0.05",
                      pCutoff = 0.05, 
                      FCcutoff = 0.58,
                      pointSize = 4,
                      labSize = 4)

cowplot::plot_grid(p3)



#### 2. SIRT4 #####
# Filter for down-regulated genes (negative log2FoldChange)
sirt4_down_genes <- unique_genes_full_list$SIRT4_KO %>% 
  filter(log2FoldChange < 0) %>% 
  pull(gene_name)

sirt4_down_g <- downregulated_sig_lists$SIRT4 %>%
  pull(gene_name)

# Convert SIRT4 down-regulated genes to Entrez IDs
sirt4_down_entrez_ids <- bitr(sirt4_down_genes,
                              fromType = "SYMBOL",
                              toType = "ENTREZID",
                              OrgDb = org.Hs.eg.db)$ENTREZID

length(sirt4_down_entrez_ids) # 29

sirt4_down_g_entrez_ids <- bitr(
  sirt4_down_g,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)$ENTREZID

# Target Entrez IDs to extract (These are genes from DNA damage response unique GO pathways)
target_entrez <- c("79913", "9212", "580", "596", "602", "641", "90427", "672", "675",
                   "993", "1032", "1293", "51514", "1894", "2138", "2177", "2189", "2237",
                   "55693", "3815", "4436", "10276", "81832", "64710", "55872", "5111", "8863",
                   "5424", "10714", "5810", "5888", "10635", "25788", "8438", "55159", "83695",
                   "6240", "51435", "23657", "7056", "90381", "54962", "11073", "7398", "7516", "7518")
# Extract only these genes from your vector
sirt4_extracted_entrez <- sirt4_down_g_entrez_ids[sirt4_down_g_entrez_ids %in% target_entrez]

# --- 2. Get DNA repair gene sets from MSigDB ---
# Get all human gene sets from Hallmark and C2 (Reactome) collections
msig_h <- msigdbr(species = "Homo sapiens", collection = "H")
msig_c2 <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")

# Check the column names to see what's available
print("MSigDB column names:")
print(colnames(msig_h)) # we need ncbi_gene column that contains ENTREZ ID 

# Filter for DNA repair-related pathways
dna_repair_terms <- c(
  # Hallmark DNA Repair
  "HALLMARK_DNA_REPAIR",
  
  # Reactome DNA Repair Pathways
  "REACTOME_HOMOLOGY_DIRECTED_REPAIR",
  "REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR",
  "REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ",
  "REACTOME_BASE_EXCISION_REPAIR",
  "REACTOME_NUCLEOTIDE_EXCISION_REPAIR",
  "REACTOME_MISMATCH_REPAIR",
  "REACTOME_GLOBAL_GENOME_NUCLEOTIDE_EXCISION_REPAIR_GG_NER",
  "REACTOME_TRANSCRIPTION_COUPLED_NUCLEOTIDE_EXCISION_REPAIR_TC_NER",
  "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR",
  "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS"
)

# Create the TERM2GENE data frame directly
term2gene <- data.frame()

for(term in dna_repair_terms) {
  if(term %in% unique(msig_h$gs_name)) {
    genes_df <- msig_h %>% 
      filter(gs_name == term) %>%
      dplyr::select(term = gs_name, gene = ncbi_gene)
  } else if(term %in% unique(msig_c2$gs_name)) {
    genes_df <- msig_c2 %>% 
      filter(gs_name == term) %>%
      dplyr::select(term = gs_name, gene = ncbi_gene)
  } else {
    warning(paste("Term not found:", term))
    next
  }
  term2gene <- rbind(term2gene, genes_df)
}

# Convert gene column to character (important for matching)
term2gene$gene <- as.character(term2gene$gene)

print("DNA repair gene sets summary:")
print(table(term2gene$term))

# --- 3. Check gene ID compatibility ---
# Check if any of our genes are in the DNA repair sets
intersection <- intersect(sirt4_extracted_entrez, term2gene$gene)
print(paste("Number of SIRT4 down genes in DNA repair pathways:", length(intersection)))

if(length(intersection) == 0) {
  print("No overlap found. Checking gene ID formats...")
  print(paste("Sample of SIRT4 Entrez IDs:", paste(head(sirt4_down_entrez_ids), collapse = ", ")))
  print(paste("Sample of DNA repair Entrez IDs:", paste(head(unique(term2gene$gene)), collapse = ", ")))
}

# --- 4. Perform enrichment analysis ---
if(length(intersection) > 0) {
  sirt4_dna_repair_enrich <- enricher(
    gene = sirt4_extracted_entrez,
    universe = background_entrez_ids,
    TERM2GENE = term2gene,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 3,    # Reduced for more sensitivity
    maxGSSize = 500
  )
  
  # --- 5. View and interpret results ---
  if(!is.null(sirt4_dna_repair_enrich)) {
    
    results_df <- sirt4_dna_repair_enrich@result
    significant_results <- results_df %>% filter(p.adjust < 0.05)
    
    print("SIRT4 Down-regulated Genes - DNA Repair Pathway Enrichment:")
    print(paste("Total pathways tested:", nrow(results_df)))
    print(paste("Significant pathways (p.adj < 0.05):", nrow(significant_results)))
    
    if(nrow(significant_results) > 0) {
      # Print significant results
      print("Significant DNA Repair Pathways:")
      print(significant_results %>% 
              dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count) %>%
              arrange(p.adjust))
      
      # Create visualization
      plot_data <- significant_results %>%
        arrange(p.adjust) %>%
        head(15)
      
      # Create dot plot - color by p.adjust directly
      p <- ggplot(plot_data, aes(x = Count, y = reorder(Description, Count), 
                                 size = Count, color = p.adjust)) +
        geom_point(alpha = 0.8) +
        scale_color_gradient(low = "red", high = "blue", 
                             name = "Adjusted p-value",
                             trans = "reverse") +  # Reverse so red = more significant
        scale_size_continuous(range = c(3, 8), name = "Gene Count") +
        labs(title = "SIRT4: DNA repair pathways enriched in down-regulated genes",
             subtitle = "Pathways normally activated by SIRT4",
             x = "Number of Genes", y = "") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5),
              axis.text.y = element_text(size = 10))
      
      print(p)
      
      # Extract and examine specific genes
      print("Genes driving significant enrichments:")
      for(i in 1:nrow(plot_data)) {
        pathway <- plot_data$ID[i]
        genes_in_pathway <- strsplit(plot_data$geneID[i], "/")[[1]]
        
        # Convert back to symbols for readability
        gene_symbols <- bitr(genes_in_pathway,
                             fromType = "ENTREZID", 
                             toType = "SYMBOL",
                             OrgDb = org.Hs.eg.db)$SYMBOL
        
        print(paste("=== ", pathway, " ==="))
        print(paste("Genes (", length(gene_symbols), "):", paste(gene_symbols, collapse = ", ")))
        cat("\n")
      }
      
    } else {
      print("No significant DNA repair pathways found.")
      print("Top pathways (even if not significant):")
      print(results_df %>% 
              arrange(p.adjust) %>% 
              head(10) %>%
              dplyr::select(ID, Description, p.adjust, Count))
    }
  }
} else {
  print("No overlap between SIRT4 down-regulated genes and DNA repair pathways.")
}



# Volcano plot for SIRT3
library(EnhancedVolcano)

# Filter to get only the 5 DNA repair genes
dna_repair_genes <- downregulated_sig_lists$SIRT4 %>%
  filter(gene_name %in% c("BARD1", "BLM", "BRCA1", "BRCA2", "EYA1", "FEN1", "PCNA", "POLD1", "POLD3", 
                          "RAD1", "RAD51", "RAD51AP1", "RHNO1", "TIPIN", "TOPBP1", "XRCC2", "XRCC4", "MSH2"))

# Simple volcano plot with just these genes
p4 <- EnhancedVolcano(dna_repair_genes,
                      lab = dna_repair_genes$gene_name, 
                      x = 'log2FoldChange',
                      y = 'padj', 
                      title = 'Down-regulated DNA repair genes in SIRT4 KO',
                      subtitle = "FC cutoff = 0.58; p-value cutoff = 0.05",
                      pCutoff = 0.05, 
                      FCcutoff = 0.58,
                      pointSize = 4,
                      labSize = 4)

cowplot::plot_grid(p4)




#### 3. SIRT5 #####
# Filter for down-regulated genes (negative log2FoldChange)
sirt5_down_genes <- unique_genes_full_list$SIRT5_KO %>% 
  filter(log2FoldChange < 0) %>% 
  pull(gene_name)

sirt5_down_g <- downregulated_sig_lists$SIRT5 %>%
  pull(gene_name)

# Convert SIRT5 down-regulated genes to Entrez IDs
sirt5_down_entrez_ids <- bitr(sirt5_down_genes,
                              fromType = "SYMBOL",
                              toType = "ENTREZID",
                              OrgDb = org.Hs.eg.db)$ENTREZID

length(sirt5_down_entrez_ids) # 20

sirt5_down_g_entrez_ids <- bitr(
  sirt5_down_g,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)$ENTREZID

# Target Entrez IDs to extract (These are genes from DNA damage response unique GO pathways)
target_entrez <- c("79000", "2187", "92797", "8091", "10459", "55010", "5347", 
                   "10721", "116028", "9656", "83695", "6117", "84250", "11073", "7518")

# Extract only these genes from your vector
sirt5_extracted_entrez <- sirt5_down_g_entrez_ids[sirt5_down_g_entrez_ids %in% target_entrez]

# --- 2. Get DNA repair gene sets from MSigDB ---
# Get all human gene sets from Hallmark and C2 (Reactome) collections
msig_h <- msigdbr(species = "Homo sapiens", collection = "H")
msig_c2 <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")

# Check the column names to see what's available
print("MSigDB column names:")
print(colnames(msig_h)) # we need ncbi_gene column that contains ENTREZ ID 

# Filter for DNA repair-related pathways
dna_repair_terms <- c(
  # Hallmark DNA Repair
  "HALLMARK_DNA_REPAIR",
  
  # Reactome DNA Repair Pathways
  "REACTOME_HOMOLOGY_DIRECTED_REPAIR",
  "REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR",
  "REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ",
  "REACTOME_BASE_EXCISION_REPAIR",
  "REACTOME_NUCLEOTIDE_EXCISION_REPAIR",
  "REACTOME_MISMATCH_REPAIR",
  "REACTOME_GLOBAL_GENOME_NUCLEOTIDE_EXCISION_REPAIR_GG_NER",
  "REACTOME_TRANSCRIPTION_COUPLED_NUCLEOTIDE_EXCISION_REPAIR_TC_NER",
  "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR",
  "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS"
)

# Create the TERM2GENE data frame directly
term2gene <- data.frame()

for(term in dna_repair_terms) {
  if(term %in% unique(msig_h$gs_name)) {
    genes_df <- msig_h %>% 
      filter(gs_name == term) %>%
      dplyr::select(term = gs_name, gene = ncbi_gene)
  } else if(term %in% unique(msig_c2$gs_name)) {
    genes_df <- msig_c2 %>% 
      filter(gs_name == term) %>%
      dplyr::select(term = gs_name, gene = ncbi_gene)
  } else {
    warning(paste("Term not found:", term))
    next
  }
  term2gene <- rbind(term2gene, genes_df)
}

# Convert gene column to character (important for matching)
term2gene$gene <- as.character(term2gene$gene)

print("DNA repair gene sets summary:")
print(table(term2gene$term))

# --- 3. Check gene ID compatibility ---
# Check if any of our genes are in the DNA repair sets
intersection <- intersect(sirt5_extracted_entrez, term2gene$gene)
print(paste("Number of SIRT5 down genes in DNA repair pathways:", length(intersection)))

if(length(intersection) == 0) {
  print("No overlap found. Checking gene ID formats...")
  print(paste("Sample of SIRT5 Entrez IDs:", paste(head(sirt5_down_entrez_ids), collapse = ", ")))
  print(paste("Sample of DNA repair Entrez IDs:", paste(head(unique(term2gene$gene)), collapse = ", ")))
}

# --- 4. Perform enrichment analysis ---
if(length(intersection) > 0) {
  sirt5_dna_repair_enrich <- enricher(
    gene = sirt5_extracted_entrez,
    universe = background_entrez_ids,
    TERM2GENE = term2gene,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 3,    # Reduced for more sensitivity
    maxGSSize = 500
  )
  
  # --- 5. View and interpret results ---
  if(!is.null(sirt5_dna_repair_enrich)) {
    
    results_df <- sirt5_dna_repair_enrich@result
    significant_results <- results_df %>% filter(p.adjust < 0.05)
    
    print("SIRT5 Down-regulated Genes - DNA Repair Pathway Enrichment:")
    print(paste("Total pathways tested:", nrow(results_df)))
    print(paste("Significant pathways (p.adj < 0.05):", nrow(significant_results)))
    
    if(nrow(significant_results) > 0) {
      # Print significant results
      print("Significant DNA Repair Pathways:")
      print(significant_results %>% 
              dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count) %>%
              arrange(p.adjust))
      
      # Create visualization
      plot_data <- significant_results %>%
        arrange(p.adjust) %>%
        head(15)
      
      # Create dot plot - color by p.adjust directly
      p <- ggplot(plot_data, aes(x = Count, y = reorder(Description, Count), 
                                 size = Count, color = p.adjust)) +
        geom_point(alpha = 0.8) +
        scale_color_gradient(low = "red", high = "blue", 
                             name = "Adjusted p-value",
                             trans = "reverse") +  # Reverse so red = more significant
        scale_size_continuous(range = c(3, 8), name = "Gene Count") +
        labs(title = "SIRT5: DNA repair pathways enriched in down-regulated genes",
             subtitle = "Pathways normally activated by SIRT5",
             x = "Number of Genes", y = "") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5),
              axis.text.y = element_text(size = 10))
      
      print(p)
      
      # Extract and examine specific genes
      print("Genes driving significant enrichments:")
      for(i in 1:nrow(plot_data)) {
        pathway <- plot_data$ID[i]
        genes_in_pathway <- strsplit(plot_data$geneID[i], "/")[[1]]
        
        # Convert back to symbols for readability
        gene_symbols <- bitr(genes_in_pathway,
                             fromType = "ENTREZID", 
                             toType = "SYMBOL",
                             OrgDb = org.Hs.eg.db)$SYMBOL
        
        print(paste("=== ", pathway, " ==="))
        print(paste("Genes (", length(gene_symbols), "):", paste(gene_symbols, collapse = ", ")))
        cat("\n")
      }
      
    } else {
      print("No significant DNA repair pathways found.")
      print("Top pathways (even if not significant):")
      print(results_df %>% 
              arrange(p.adjust) %>% 
              head(10) %>%
              dplyr::select(ID, Description, p.adjust, Count))
    }
  }
} else {
  print("No overlap between SIRT5 down-regulated genes and DNA repair pathways.")
}


# 79728 - SIRT5 gene that is intersected with DNA repair pathways 
gene_symbol <- bitr("79728",
                    fromType = "ENTREZID",
                    toType = "SYMBOL",
                    OrgDb = org.Hs.eg.db)

# Volcano plot 
library(EnhancedVolcano)
# Create the keyvals color vector with names
keyvals <- ifelse(
  unique_genes_full_list$SIRT5$gene_name == 'PALB2', 'red', 'grey60'
)
names(keyvals)[keyvals == 'red'] <- 'PALB2'
names(keyvals)[keyvals == 'grey60'] <- 'other_genes'

# Create the volcano plot
p5 <- EnhancedVolcano(unique_genes_full_list$SIRT5,
                      lab = unique_genes_full_list$SIRT5$gene_name,
                      x = "log2FoldChange",
                      y = "padj",
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      pointSize = ifelse(unique_genes_full_list$SIRT5$gene_name == 'PALB2', 4, 2),
                      labSize = 3.5,
                      selectLab = 'PALB2',
                      title = "Unique SIRT5_KO genes vs WT",
                      subtitle = "FC cutoff = 0.58; p-value cutoff = 0.05",
                      colCustom = keyvals,
                      colAlpha = 0.8,
                      drawConnectors = TRUE,
                      widthConnectors = 0.2,
                      legendPosition = 'right')

cowplot::plot_grid(p5)



# Volcano plot for SIRT3
# Filter to get only the 5 DNA repair genes
dna_repair_genes <- downregulated_sig_lists$SIRT5 %>%
  filter(gene_name %in% c('MDC1', 'POLQ', 'RHNO1', 'RMI2', 'RPA1', 'TOPBP1', 'XRCC4'))

# Simple volcano plot with just these genes
p5 <- EnhancedVolcano(dna_repair_genes,
                      lab = dna_repair_genes$gene_name, 
                      x = 'log2FoldChange',
                      y = 'padj', 
                      title = 'Down-regulated DNA repair genes in SIRT5 KO',
                      subtitle = "FC cutoff = 0.58; p-value cutoff = 0.05",
                      pCutoff = 0.05, 
                      FCcutoff = 0.58,
                      pointSize = 4,
                      labSize = 4)

cowplot::plot_grid(p5)




#### VII. Reactome enrichment analysis using compareCluster function ####

##### a. For up-regulated genes ######

# 1. Run the comparative analysis for Reactome pathways
# Filter out empty lists to prevent errors
upregulated_to_compare <- upregulated_entrez_lists[sapply(upregulated_entrez_lists, length) > 0]

compare_reactome_results <- compareCluster(
  geneCluster   = upregulated_to_compare,
  fun           = "enrichPathway",      # <-- Use enrichPathway for Reactome
  universe      = background_entrez_ids,
  organism      = "human",              # Argument for enrichPathway
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1
)

# 2. Filter and visualize only unique pathways 
# First, check if the analysis returned any results
if (!is.null(compare_reactome_results) && nrow(compare_reactome_results) > 0) {
  
  # 1. Convert the result object to a data.frame for easy filtering
  results_df_reactome <- as.data.frame(compare_reactome_results)
  
  # 2. Find the pathways that are unique to a single cluster (SIRT)
  unique_pathways_df_reactome <- results_df_reactome %>%
    group_by(Description) %>%
    filter(n() == 1) %>%
    ungroup()
  
  # 3. Create a new compareCluster object containing ONLY the unique pathways
  compare_reactome_unique <- compare_reactome_results
  compare_reactome_unique@compareClusterResult <- unique_pathways_df_reactome
  
  # 4. Check if any unique pathways were found before proceeding
  if(nrow(compare_reactome_unique) > 0) {
    
    # --- 5. Modify the cluster names from "SIRT1_KO" to "1", etc. ---
    results_to_plot_df <- as.data.frame(compare_reactome_unique)
    results_to_plot_df$Cluster <- gsub("SIRT|_KO", "", results_to_plot_df$Cluster)
    results_to_plot_df$Cluster <- factor(results_to_plot_df$Cluster, levels = as.character(1:7))
    
    # Put the modified data back into the object for plotting
    compare_reactome_unique@compareClusterResult <- results_to_plot_df
    
    # --- 6. Create the final dot plot ---
    p_unique_reactome <- dotplot(
      compare_reactome_unique, 
      showCategory = 2,      # <-- As requested, showing top 2 unique for each
      by = "p.adjust"
    ) +
      ggtitle("Unique Reactome enrichment across up- genes") +
      xlab("SIRT KO") +
      ylab("")
    
    print(p_unique_reactome)
    
  } else {
    cat("After filtering, no unique Reactome pathways were found to display.\n")
  }
  
} else {
  cat("The initial compareCluster analysis for Reactome returned no significant results.\n")
}

# Split the unique pathways table into a list of tables
list_of_unique_reactome_tables <- split(unique_pathways_df_reactome, unique_pathways_df_reactome$Cluster)




######## For SIRT1 
sirt1_unique_reactome_df <- list_of_unique_reactome_tables$SIRT1_KO

if (!is.null(sirt1_unique_reactome_df) && nrow(sirt1_unique_reactome_df) > 0) {
  
  # --- Export the table for manual thematic grouping ---
  write.csv(sirt1_unique_reactome_df, "SIRT1_unique_up_Reactome_for_grouping.csv", row.names = FALSE)
  
  cat("Table 'SIRT1_unique_up_Reactome_for_grouping.csv' was saved.\n")
  cat("Please open it in a spreadsheet editor, add a 'Theme' column, and save it as a new file.\n")
  
} else {
  cat("No unique Reactome pathways found for SIRT1_KO.\n")
}

# Manually creat the "grouped" CSV file 

# Read your new, grouped table 
grouped_data_reactome <- read.csv2("SIRT1_unique_up_Reactome_grouped.csv", sep=";")

# Prepare the data for plotting 
if (nrow(grouped_data_reactome) > 0) {
  
  plot_data_reactome <- grouped_data_reactome %>%
    mutate(
      p.adjust = as.numeric(gsub(",", ".", p.adjust)),
      Count = as.numeric(Count)
    ) %>%
    filter(!is.na(p.adjust), !is.na(Theme))
  
  plot_data_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice(1:5) %>% # Show top 5 pathways per theme
    ungroup()
  
  theme_order_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    summarise(max_count = max(Count)) %>%
    arrange(desc(max_count)) %>%
    pull(Theme)
  
  plot_data_reactome <- plot_data_reactome %>%
    mutate(
      Description = reorder_within(Description, Count, Theme),
      Theme = factor(Theme, levels = rev(theme_order_reactome))
    )
  
  # --- Create the customized bar plot ---
  p_final_reactome <- ggplot(plot_data_reactome, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_y_reordered() +
    scale_fill_gradientn(colors = c("red", "magenta", "blue"), name = "p.adjust") +
    labs(title = "Unique Reactome Pathways for SIRT1_KO Upregulated Genes", x = "Count", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
      strip.background = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  # --- Display the final plot ---
  print(p_final_reactome)
}


# Extract the genes contributing to unique Reactome pathways 
# 1. Get the names of the pathways that are unique to SIRT1_KO
sirt1_unique_reactome_names <- unique_pathways_df_reactome %>%
  dplyr::filter(Cluster == "SIRT1_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df_reactome <- as.data.frame(compare_reactome_results)

# 3. Filter the full results
sirt1_unique_genes_df_reactome <- full_results_df_reactome %>%
  dplyr::filter(
    Cluster == "SIRT1_KO",
    Description %in% sirt1_unique_reactome_names
  )

# 4. Extract, split, and combine all the gene IDs
all_sirt1_unique_genes_reactome <- sirt1_unique_genes_df_reactome %>%
  pull(geneID) %>%
  strsplit(split = "/") %>%
  unlist() %>%
  unique()

# 5. Count them!
number_of_sirt1_genes_reactome <- length(all_sirt1_unique_genes_reactome)

cat("The number of genes contributing to the unique Reactome functions of SIRT1 is:",
    number_of_sirt1_genes_reactome, "\n")



######## For SIRT2
sirt2_unique_reactome_df <- list_of_unique_reactome_tables$SIRT2_KO

if (!is.null(sirt2_unique_reactome_df) && nrow(sirt2_unique_reactome_df) > 0) {
  
  # --- Export the table for manual thematic grouping ---
  write.csv(sirt2_unique_reactome_df, "SIRT2_unique_up_Reactome_for_grouping.csv", row.names = FALSE)
  
  cat("Table 'SIRT2_unique_up_Reactome_for_grouping.csv' was saved.\n")
  cat("Please open it in a spreadsheet editor, add a 'Theme' column, and save it as a new file.\n")
  
} else {
  cat("No unique Reactome pathways found for SIRT2_KO.\n")
}

# Manually creat the "grouped" CSV file 

# Read your new, grouped table 
grouped_data_reactome <- read.csv2("SIRT2_unique_up_Reactome_grouped.csv", sep=";")

# Prepare the data for plotting 
if (nrow(grouped_data_reactome) > 0) {
  
  plot_data_reactome <- grouped_data_reactome %>%
    mutate(
      p.adjust = as.numeric(gsub(",", ".", p.adjust)),
      Count = as.numeric(Count)
    ) %>%
    filter(!is.na(p.adjust), !is.na(Theme))
  
  plot_data_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice(1:5) %>% # Show top 5 pathways per theme
    ungroup()
  
  theme_order_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    summarise(max_count = max(Count)) %>%
    arrange(desc(max_count)) %>%
    pull(Theme)
  
  plot_data_reactome <- plot_data_reactome %>%
    mutate(
      Description = reorder_within(Description, Count, Theme),
      Theme = factor(Theme, levels = rev(theme_order_reactome))
    )
  
  # --- Create the customized bar plot ---
  p_final_reactome <- ggplot(plot_data_reactome, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_y_reordered() +
    scale_fill_gradientn(colors = c("red", "magenta", "blue"), name = "p.adjust") +
    labs(title = "Unique Reactome Pathways for SIRT2_KO Upregulated Genes", x = "Count", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
      strip.background = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  # --- Display the final plot ---
  print(p_final_reactome)
}


# Extract the genes contributing to unique Reactome pathways 
# 1. Get the names of the pathways that are unique to SIRT2_KO
sirt2_unique_reactome_names <- unique_pathways_df_reactome %>%
  dplyr::filter(Cluster == "SIRT2_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df_reactome <- as.data.frame(compare_reactome_results)

# 3. Filter the full results
sirt2_unique_genes_df_reactome <- full_results_df_reactome %>%
  dplyr::filter(
    Cluster == "SIRT2_KO",
    Description %in% sirt2_unique_reactome_names
  )

# 4. Extract, split, and combine all the gene IDs
all_sirt2_unique_genes_reactome <- sirt2_unique_genes_df_reactome %>%
  pull(geneID) %>%
  strsplit(split = "/") %>%
  unlist() %>%
  unique()

# 5. Count them!
number_of_sirt2_genes_reactome <- length(all_sirt2_unique_genes_reactome)

cat("The number of genes contributing to the unique Reactome functions of SIRT2 is:",
    number_of_sirt2_genes_reactome, "\n")



######## For SIRT3
sirt3_unique_reactome_df <- list_of_unique_reactome_tables$SIRT3_KO

if (!is.null(sirt3_unique_reactome_df) && nrow(sirt3_unique_reactome_df) > 0) {
  
  # --- Export the table for manual thematic grouping ---
  write.csv(sirt3_unique_reactome_df, "SIRT3_unique_up_Reactome_for_grouping.csv", row.names = FALSE)
  
  cat("Table 'SIRT3_unique_up_Reactome_for_grouping.csv' was saved.\n")
  cat("Please open it in a spreadsheet editor, add a 'Theme' column, and save it as a new file.\n")
  
} else {
  cat("No unique Reactome pathways found for SIRT3_KO.\n")
}

# Manually creat the "grouped" CSV file 

# Read your new, grouped table 
grouped_data_reactome <- read.csv2("SIRT3_unique_up_Reactome_grouped.csv", sep=";")

# Prepare the data for plotting 
if (nrow(grouped_data_reactome) > 0) {
  
  plot_data_reactome <- grouped_data_reactome %>%
    mutate(
      p.adjust = as.numeric(gsub(",", ".", p.adjust)),
      Count = as.numeric(Count)
    ) %>%
    filter(!is.na(p.adjust), !is.na(Theme))
  
  plot_data_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice(1:5) %>% # Show top 5 pathways per theme
    ungroup()
  
  theme_order_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    summarise(max_count = max(Count)) %>%
    arrange(desc(max_count)) %>%
    pull(Theme)
  
  plot_data_reactome <- plot_data_reactome %>%
    mutate(
      Description = reorder_within(Description, Count, Theme),
      Theme = factor(Theme, levels = rev(theme_order_reactome))
    )
  
  # --- Create the customized bar plot ---
  p_final_reactome <- ggplot(plot_data_reactome, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_y_reordered() +
    scale_fill_gradientn(colors = c("red", "magenta", "blue"), name = "p.adjust") +
    labs(title = "Unique Reactome Pathways for SIRT3_KO Upregulated Genes", x = "Count", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
      strip.background = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  # --- Display the final plot ---
  print(p_final_reactome)
}


# Extract the genes contributing to unique Reactome pathways 
# 1. Get the names of the pathways that are unique to SIRT3_KO
sirt3_unique_reactome_names <- unique_pathways_df_reactome %>%
  dplyr::filter(Cluster == "SIRT3_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df_reactome <- as.data.frame(compare_reactome_results)

# 3. Filter the full results
sirt3_unique_genes_df_reactome <- full_results_df_reactome %>%
  dplyr::filter(
    Cluster == "SIRT3_KO",
    Description %in% sirt3_unique_reactome_names
  )

# 4. Extract, split, and combine all the gene IDs
all_sirt3_unique_genes_reactome <- sirt3_unique_genes_df_reactome %>%
  pull(geneID) %>%
  strsplit(split = "/") %>%
  unlist() %>%
  unique()

# 5. Count them!
number_of_sirt3_genes_reactome <- length(all_sirt3_unique_genes_reactome)

cat("The number of genes contributing to the unique Reactome functions of SIRT3 is:",
    number_of_sirt3_genes_reactome, "\n")




######## For SIRT4
sirt4_unique_reactome_df <- list_of_unique_reactome_tables$SIRT4_KO

if (!is.null(sirt4_unique_reactome_df) && nrow(sirt4_unique_reactome_df) > 0) {
  
  # --- Export the table for manual thematic grouping ---
  write.csv(sirt4_unique_reactome_df, "SIRT4_unique_up_Reactome_for_grouping.csv", row.names = FALSE)
  
  cat("Table 'SIRT4_unique_up_Reactome_for_grouping.csv' was saved.\n")
  cat("Please open it in a spreadsheet editor, add a 'Theme' column, and save it as a new file.\n")
  
} else {
  cat("No unique Reactome pathways found for SIRT4_KO.\n")
}

# !No unique Reactome pathways found for SIRT4_KO!




######## For SIRT5
sirt5_unique_reactome_df <- list_of_unique_reactome_tables$SIRT5_KO

if (!is.null(sirt5_unique_reactome_df) && nrow(sirt5_unique_reactome_df) > 0) {
  
  # --- Export the table for manual thematic grouping ---
  write.csv(sirt5_unique_reactome_df, "SIRT5_unique_up_Reactome_for_grouping.csv", row.names = FALSE)
  
  cat("Table 'SIRT5_unique_up_Reactome_for_grouping.csv' was saved.\n")
  cat("Please open it in a spreadsheet editor, add a 'Theme' column, and save it as a new file.\n")
  
} else {
  cat("No unique Reactome pathways found for SIRT5_KO.\n")
}

# Manually creat the "grouped" CSV file 

# Read your new, grouped table 
grouped_data_reactome <- read.csv2("SIRT5_unique_up_Reactome_grouped.csv", sep=";")

# Prepare the data for plotting 
if (nrow(grouped_data_reactome) > 0) {
  
  plot_data_reactome <- grouped_data_reactome %>%
    mutate(
      p.adjust = as.numeric(gsub(",", ".", p.adjust)),
      Count = as.numeric(Count)
    ) %>%
    filter(!is.na(p.adjust), !is.na(Theme))
  
  plot_data_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice(1:5) %>% # Show top 5 pathways per theme
    ungroup()
  
  theme_order_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    summarise(max_count = max(Count)) %>%
    arrange(desc(max_count)) %>%
    pull(Theme)
  
  plot_data_reactome <- plot_data_reactome %>%
    mutate(
      Description = reorder_within(Description, Count, Theme),
      Theme = factor(Theme, levels = rev(theme_order_reactome))
    )
  
  # --- Create the customized bar plot ---
  p_final_reactome <- ggplot(plot_data_reactome, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_y_reordered() +
    scale_fill_gradientn(colors = c("red", "magenta", "blue"), name = "p.adjust") +
    labs(title = "Unique Reactome Pathways for SIRT5_KO Upregulated Genes", x = "Count", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
      strip.background = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  # --- Display the final plot ---
  print(p_final_reactome)
}


# Extract the genes contributing to unique Reactome pathways 
# 1. Get the names of the pathways that are unique to SIRT5_KO
sirt5_unique_reactome_names <- unique_pathways_df_reactome %>%
  dplyr::filter(Cluster == "SIRT5_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df_reactome <- as.data.frame(compare_reactome_results)

# 3. Filter the full results
sirt5_unique_genes_df_reactome <- full_results_df_reactome %>%
  dplyr::filter(
    Cluster == "SIRT5_KO",
    Description %in% sirt5_unique_reactome_names
  )

# 4. Extract, split, and combine all the gene IDs
all_sirt5_unique_genes_reactome <- sirt5_unique_genes_df_reactome %>%
  pull(geneID) %>%
  strsplit(split = "/") %>%
  unlist() %>%
  unique()

# 5. Count them!
number_of_sirt5_genes_reactome <- length(all_sirt5_unique_genes_reactome)

cat("The number of genes contributing to the unique Reactome functions of SIRT5 is:",
    number_of_sirt5_genes_reactome, "\n")



######## For SIRT6
sirt6_unique_reactome_df <- list_of_unique_reactome_tables$SIRT6_KO

if (!is.null(sirt6_unique_reactome_df) && nrow(sirt6_unique_reactome_df) > 0) {
  
  # --- Export the table for manual thematic grouping ---
  write.csv(sirt6_unique_reactome_df, "SIRT6_unique_up_Reactome_for_grouping.csv", row.names = FALSE)
  
  cat("Table 'SIRT6_unique_up_Reactome_for_grouping.csv' was saved.\n")
  cat("Please open it in a spreadsheet editor, add a 'Theme' column, and save it as a new file.\n")
  
} else {
  cat("No unique Reactome pathways found for SIRT6_KO.\n")
}

# Manually creat the "grouped" CSV file 

# Read your new, grouped table 
grouped_data_reactome <- read.csv2("SIRT6_unique_up_Reactome_grouped.csv", sep=";")

# Prepare the data for plotting 
if (nrow(grouped_data_reactome) > 0) {
  
  plot_data_reactome <- grouped_data_reactome %>%
    mutate(
      p.adjust = as.numeric(gsub(",", ".", p.adjust)),
      Count = as.numeric(Count)
    ) %>%
    filter(!is.na(p.adjust), !is.na(Theme))
  
  plot_data_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice(1:5) %>% # Show top 5 pathways per theme
    ungroup()
  
  theme_order_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    summarise(max_count = max(Count)) %>%
    arrange(desc(max_count)) %>%
    pull(Theme)
  
  plot_data_reactome <- plot_data_reactome %>%
    mutate(
      Description = reorder_within(Description, Count, Theme),
      Theme = factor(Theme, levels = rev(theme_order_reactome))
    )
  
  # --- Create the customized bar plot ---
  p_final_reactome <- ggplot(plot_data_reactome, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_y_reordered() +
    scale_fill_gradientn(colors = c("red", "magenta", "blue"), name = "p.adjust") +
    labs(title = "Unique Reactome Pathways for SIRT6_KO Upregulated Genes", x = "Count", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
      strip.background = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  # --- Display the final plot ---
  print(p_final_reactome)
}


# Extract the genes contributing to unique Reactome pathways 
# 1. Get the names of the pathways that are unique to SIRT6_KO
sirt6_unique_reactome_names <- unique_pathways_df_reactome %>%
  dplyr::filter(Cluster == "SIRT6_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df_reactome <- as.data.frame(compare_reactome_results)

# 3. Filter the full results
sirt6_unique_genes_df_reactome <- full_results_df_reactome %>%
  dplyr::filter(
    Cluster == "SIRT6_KO",
    Description %in% sirt6_unique_reactome_names
  )

# 4. Extract, split, and combine all the gene IDs
all_sirt6_unique_genes_reactome <- sirt6_unique_genes_df_reactome %>%
  pull(geneID) %>%
  strsplit(split = "/") %>%
  unlist() %>%
  unique()

# 5. Count them!
number_of_sirt6_genes_reactome <- length(all_sirt6_unique_genes_reactome)

cat("The number of genes contributing to the unique Reactome functions of SIRT6 is:",
    number_of_sirt6_genes_reactome, "\n")




######## For SIRT7
sirt7_unique_reactome_df <- list_of_unique_reactome_tables$SIRT7_KO

if (!is.null(sirt7_unique_reactome_df) && nrow(sirt7_unique_reactome_df) > 0) {
  
  # --- Export the table for manual thematic grouping ---
  write.csv(sirt7_unique_reactome_df, "SIRT7_unique_up_Reactome_for_grouping.csv", row.names = FALSE)
  
  cat("Table 'SIRT7_unique_up_Reactome_for_grouping.csv' was saved.\n")
  cat("Please open it in a spreadsheet editor, add a 'Theme' column, and save it as a new file.\n")
  
} else {
  cat("No unique Reactome pathways found for SIRT7_KO.\n")
}

# Manually creat the "grouped" CSV file 

# Read your new, grouped table 
grouped_data_reactome <- read.csv2("SIRT7_unique_up_Reactome_grouped.csv", sep=";")

# Prepare the data for plotting 
if (nrow(grouped_data_reactome) > 0) {
  
  plot_data_reactome <- grouped_data_reactome %>%
    mutate(
      p.adjust = as.numeric(gsub(",", ".", p.adjust)),
      Count = as.numeric(Count)
    ) %>%
    filter(!is.na(p.adjust), !is.na(Theme))
  
  plot_data_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice(1:5) %>% # Show top 5 pathways per theme
    ungroup()
  
  theme_order_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    summarise(max_count = max(Count)) %>%
    arrange(desc(max_count)) %>%
    pull(Theme)
  
  plot_data_reactome <- plot_data_reactome %>%
    mutate(
      Description = reorder_within(Description, Count, Theme),
      Theme = factor(Theme, levels = rev(theme_order_reactome))
    )
  
  # --- Create the customized bar plot ---
  p_final_reactome <- ggplot(plot_data_reactome, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_y_reordered() +
    scale_fill_gradientn(colors = c("red", "magenta", "blue"), name = "p.adjust") +
    labs(title = "Unique Reactome Pathways for SIRT7_KO Upregulated Genes", x = "Count", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
      strip.background = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  # --- Display the final plot ---
  print(p_final_reactome)
}


# Extract the genes contributing to unique Reactome pathways 
# 1. Get the names of the pathways that are unique to SIRT7_KO
sirt7_unique_reactome_names <- unique_pathways_df_reactome %>%
  dplyr::filter(Cluster == "SIRT7_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df_reactome <- as.data.frame(compare_reactome_results)

# 3. Filter the full results
sirt7_unique_genes_df_reactome <- full_results_df_reactome %>%
  dplyr::filter(
    Cluster == "SIRT7_KO",
    Description %in% sirt7_unique_reactome_names
  )

# 4. Extract, split, and combine all the gene IDs
all_sirt7_unique_genes_reactome <- sirt7_unique_genes_df_reactome %>%
  pull(geneID) %>%
  strsplit(split = "/") %>%
  unlist() %>%
  unique()

# 5. Count them!
number_of_sirt7_genes_reactome <- length(all_sirt7_unique_genes_reactome)

cat("The number of genes contributing to the unique Reactome functions of SIRT7 is:",
    number_of_sirt7_genes_reactome, "\n")




##### b. For down-regulated genes ######

# 1. Run the comparative analysis for Reactome pathways
# Filter out empty lists to prevent errors
downregulated_to_compare <- downregulated_entrez_lists[sapply(downregulated_entrez_lists, length) > 0]

compare_reactome_results <- compareCluster(
  geneCluster   = downregulated_to_compare,
  fun           = "enrichPathway",      # <-- Use enrichPathway for Reactome
  universe      = background_entrez_ids,
  organism      = "human",              # Argument for enrichPathway
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1
)

# 2. Filter and visualize only unique pathways 
# First, check if the analysis returned any results
if (!is.null(compare_reactome_results) && nrow(compare_reactome_results) > 0) {
  
  # 1. Convert the result object to a data.frame for easy filtering
  results_df_reactome <- as.data.frame(compare_reactome_results)
  
  # 2. Find the pathways that are unique to a single cluster (SIRT)
  unique_pathways_df_reactome <- results_df_reactome %>%
    group_by(Description) %>%
    filter(n() == 1) %>%
    ungroup()
  
  # 3. Create a new compareCluster object containing ONLY the unique pathways
  compare_reactome_unique <- compare_reactome_results
  compare_reactome_unique@compareClusterResult <- unique_pathways_df_reactome
  
  # 4. Check if any unique pathways were found before proceeding
  if(nrow(compare_reactome_unique) > 0) {
    
    # --- 5. Modify the cluster names from "SIRT1_KO" to "1", etc. ---
    results_to_plot_df <- as.data.frame(compare_reactome_unique)
    results_to_plot_df$Cluster <- gsub("SIRT|_KO", "", results_to_plot_df$Cluster)
    results_to_plot_df$Cluster <- factor(results_to_plot_df$Cluster, levels = as.character(1:7))
    
    # Put the modified data back into the object for plotting
    compare_reactome_unique@compareClusterResult <- results_to_plot_df
    
    # --- 6. Create the final dot plot ---
    p_unique_reactome <- dotplot(
      compare_reactome_unique, 
      showCategory = 2,      # <-- As requested, showing top 2 unique for each
      by = "p.adjust"
    ) +
      ggtitle("Unique Reactome enrichment across down- genes") +
      xlab("SIRT KO") +
      ylab("")
    
    print(p_unique_reactome)
    
  } else {
    cat("After filtering, no unique Reactome pathways were found to display.\n")
  }
  
} else {
  cat("The initial compareCluster analysis for Reactome returned no significant results.\n")
}

# Split the unique pathways table into a list of tables
list_of_unique_reactome_tables <- split(unique_pathways_df_reactome, unique_pathways_df_reactome$Cluster)




######## For SIRT1 
sirt1_unique_reactome_df <- list_of_unique_reactome_tables$SIRT1_KO

if (!is.null(sirt1_unique_reactome_df) && nrow(sirt1_unique_reactome_df) > 0) {
  
  # --- Export the table for manual thematic grouping ---
  write.csv(sirt1_unique_reactome_df, "SIRT1_unique_down_Reactome_for_grouping.csv", row.names = FALSE)
  
  cat("Table 'SIRT1_unique_down_Reactome_for_grouping.csv' was saved.\n")
  cat("Please open it in a spreadsheet editor, add a 'Theme' column, and save it as a new file.\n")
  
} else {
  cat("No unique Reactome pathways found for SIRT1_KO.\n")
}

# Manually creat the "grouped" CSV file 

# Read your new, grouped table 
grouped_data_reactome <- read.csv2("SIRT1_unique_down_Reactome_grouped.csv", sep=";")

# Prepare the data for plotting 
if (nrow(grouped_data_reactome) > 0) {
  
  plot_data_reactome <- grouped_data_reactome %>%
    mutate(
      p.adjust = as.numeric(gsub(",", ".", p.adjust)),
      Count = as.numeric(Count)
    ) %>%
    filter(!is.na(p.adjust), !is.na(Theme))
  
  plot_data_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice(1:5) %>% # Show top 5 pathways per theme
    ungroup()
  
  theme_order_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    summarise(max_count = max(Count)) %>%
    arrange(desc(max_count)) %>%
    pull(Theme)
  
  plot_data_reactome <- plot_data_reactome %>%
    mutate(
      Description = reorder_within(Description, Count, Theme),
      Theme = factor(Theme, levels = rev(theme_order_reactome))
    )
  
  # --- Create the customized bar plot ---
  p_final_reactome <- ggplot(plot_data_reactome, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_y_reordered() +
    scale_fill_gradientn(colors = c("red", "magenta", "blue"), name = "p.adjust") +
    labs(title = "Unique Reactome Pathways for SIRT1_KO Downregulated Genes", x = "Count", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
      strip.background = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  # --- Display the final plot ---
  print(p_final_reactome)
}


# Extract the genes contributing to unique Reactome pathways 
# 1. Get the names of the pathways that are unique to SIRT1_KO
sirt1_unique_reactome_names <- unique_pathways_df_reactome %>%
  dplyr::filter(Cluster == "SIRT1_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df_reactome <- as.data.frame(compare_reactome_results)

# 3. Filter the full results
sirt1_unique_genes_df_reactome <- full_results_df_reactome %>%
  dplyr::filter(
    Cluster == "SIRT1_KO",
    Description %in% sirt1_unique_reactome_names
  )

# 4. Extract, split, and combine all the gene IDs
all_sirt1_unique_genes_reactome <- sirt1_unique_genes_df_reactome %>%
  pull(geneID) %>%
  strsplit(split = "/") %>%
  unlist() %>%
  unique()

# 5. Count them!
number_of_sirt1_genes_reactome <- length(all_sirt1_unique_genes_reactome)

cat("The number of genes contributing to the unique Reactome functions of SIRT1 is:",
    number_of_sirt1_genes_reactome, "\n")




######## For SIRT2
sirt2_unique_reactome_df <- list_of_unique_reactome_tables$SIRT2_KO

if (!is.null(sirt2_unique_reactome_df) && nrow(sirt2_unique_reactome_df) > 0) {
  
  # --- Export the table for manual thematic grouping ---
  write.csv(sirt2_unique_reactome_df, "SIRT2_unique_down_Reactome_for_grouping.csv", row.names = FALSE)
  
  cat("Table 'SIRT2_unique_down_Reactome_for_grouping.csv' was saved.\n")
  cat("Please open it in a spreadsheet editor, add a 'Theme' column, and save it as a new file.\n")
  
} else {
  cat("No unique Reactome pathways found for SIRT2_KO.\n")
}

# Manually creat the "grouped" CSV file 

# Read your new, grouped table 
grouped_data_reactome <- read.csv2("SIRT2_unique_down_Reactome_grouped.csv", sep=";")

# Prepare the data for plotting 
if (nrow(grouped_data_reactome) > 0) {
  
  plot_data_reactome <- grouped_data_reactome %>%
    mutate(
      p.adjust = as.numeric(gsub(",", ".", p.adjust)),
      Count = as.numeric(Count)
    ) %>%
    filter(!is.na(p.adjust), !is.na(Theme))
  
  plot_data_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice(1:5) %>% # Show top 5 pathways per theme
    ungroup()
  
  theme_order_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    summarise(max_count = max(Count)) %>%
    arrange(desc(max_count)) %>%
    pull(Theme)
  
  plot_data_reactome <- plot_data_reactome %>%
    mutate(
      Description = reorder_within(Description, Count, Theme),
      Theme = factor(Theme, levels = rev(theme_order_reactome))
    )
  
  # --- Create the customized bar plot ---
  p_final_reactome <- ggplot(plot_data_reactome, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_y_reordered() +
    scale_fill_gradientn(colors = c("red", "magenta", "blue"), name = "p.adjust") +
    labs(title = "Unique Reactome Pathways for SIRT2_KO Downregulated Genes", x = "Count", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
      strip.background = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  # --- Display the final plot ---
  print(p_final_reactome)
}


# Extract the genes contributing to unique Reactome pathways 
# 1. Get the names of the pathways that are unique to SIRT2_KO
sirt2_unique_reactome_names <- unique_pathways_df_reactome %>%
  dplyr::filter(Cluster == "SIRT2_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df_reactome <- as.data.frame(compare_reactome_results)

# 3. Filter the full results
sirt2_unique_genes_df_reactome <- full_results_df_reactome %>%
  dplyr::filter(
    Cluster == "SIRT2_KO",
    Description %in% sirt2_unique_reactome_names
  )

# 4. Extract, split, and combine all the gene IDs
all_sirt2_unique_genes_reactome <- sirt2_unique_genes_df_reactome %>%
  pull(geneID) %>%
  strsplit(split = "/") %>%
  unlist() %>%
  unique()

# 5. Count them!
number_of_sirt2_genes_reactome <- length(all_sirt2_unique_genes_reactome)

cat("The number of genes contributing to the unique Reactome functions of SIRT2 is:",
    number_of_sirt2_genes_reactome, "\n")




######## For SIRT3
sirt3_unique_reactome_df <- list_of_unique_reactome_tables$SIRT3_KO

if (!is.null(sirt3_unique_reactome_df) && nrow(sirt3_unique_reactome_df) > 0) {
  
  # --- Export the table for manual thematic grouping ---
  write.csv(sirt3_unique_reactome_df, "SIRT3_unique_down_Reactome_for_grouping.csv", row.names = FALSE)
  
  cat("Table 'SIRT3_unique_down_Reactome_for_grouping.csv' was saved.\n")
  cat("Please open it in a spreadsheet editor, add a 'Theme' column, and save it as a new file.\n")
  
} else {
  cat("No unique Reactome pathways found for SIRT3_KO.\n")
}

# Manually creat the "grouped" CSV file 

# Read your new, grouped table 
grouped_data_reactome <- read.csv2("SIRT3_unique_down_Reactome_grouped.csv", sep=";")

# Prepare the data for plotting 
if (nrow(grouped_data_reactome) > 0) {
  
  plot_data_reactome <- grouped_data_reactome %>%
    mutate(
      p.adjust = as.numeric(gsub(",", ".", p.adjust)),
      Count = as.numeric(Count)
    ) %>%
    filter(!is.na(p.adjust), !is.na(Theme))
  
  plot_data_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice(1:5) %>% # Show top 5 pathways per theme
    ungroup()
  
  theme_order_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    summarise(max_count = max(Count)) %>%
    arrange(desc(max_count)) %>%
    pull(Theme)
  
  plot_data_reactome <- plot_data_reactome %>%
    mutate(
      Description = reorder_within(Description, Count, Theme),
      Theme = factor(Theme, levels = rev(theme_order_reactome))
    )
  
  # --- Create the customized bar plot ---
  p_final_reactome <- ggplot(plot_data_reactome, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_y_reordered() +
    scale_fill_gradientn(colors = c("red", "magenta", "blue"), name = "p.adjust") +
    labs(title = "Unique Reactome Pathways for SIRT3_KO Downregulated Genes", x = "Count", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
      strip.background = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  # --- Display the final plot ---
  print(p_final_reactome)
}


# Extract the genes contributing to unique Reactome pathways 
# 1. Get the names of the pathways that are unique to SIRT3_KO
sirt3_unique_reactome_names <- unique_pathways_df_reactome %>%
  dplyr::filter(Cluster == "SIRT3_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df_reactome <- as.data.frame(compare_reactome_results)

# 3. Filter the full results
sirt3_unique_genes_df_reactome <- full_results_df_reactome %>%
  dplyr::filter(
    Cluster == "SIRT3_KO",
    Description %in% sirt3_unique_reactome_names
  )

# 4. Extract, split, and combine all the gene IDs
all_sirt3_unique_genes_reactome <- sirt3_unique_genes_df_reactome %>%
  pull(geneID) %>%
  strsplit(split = "/") %>%
  unlist() %>%
  unique()

# 5. Count them!
number_of_sirt3_genes_reactome <- length(all_sirt3_unique_genes_reactome)

cat("The number of genes contributing to the unique Reactome functions of SIRT3 is:",
    number_of_sirt3_genes_reactome, "\n")



######## For SIRT4
sirt4_unique_reactome_df <- list_of_unique_reactome_tables$SIRT4_KO

if (!is.null(sirt4_unique_reactome_df) && nrow(sirt4_unique_reactome_df) > 0) {
  
  # --- Export the table for manual thematic grouping ---
  write.csv(sirt4_unique_reactome_df, "SIRT4_unique_down_Reactome_for_grouping.csv", row.names = FALSE)
  
  cat("Table 'SIRT4_unique_down_Reactome_for_grouping.csv' was saved.\n")
  cat("Please open it in a spreadsheet editor, add a 'Theme' column, and save it as a new file.\n")
  
} else {
  cat("No unique Reactome pathways found for SIRT4_KO.\n")
}

# ! No unique Reactome pathways found for SIRT4_KO !




######## For SIRT5
sirt5_unique_reactome_df <- list_of_unique_reactome_tables$SIRT5_KO

if (!is.null(sirt5_unique_reactome_df) && nrow(sirt5_unique_reactome_df) > 0) {
  
  # --- Export the table for manual thematic grouping ---
  write.csv(sirt5_unique_reactome_df, "SIRT5_unique_down_Reactome_for_grouping.csv", row.names = FALSE)
  
  cat("Table 'SIRT5_unique_down_Reactome_for_grouping.csv' was saved.\n")
  cat("Please open it in a spreadsheet editor, add a 'Theme' column, and save it as a new file.\n")
  
} else {
  cat("No unique Reactome pathways found for SIRT5_KO.\n")
}

# Manually creat the "grouped" CSV file 

# Read your new, grouped table 
grouped_data_reactome <- read.csv2("SIRT5_unique_down_Reactome_grouped.csv", sep=";")

# Prepare the data for plotting 
if (nrow(grouped_data_reactome) > 0) {
  
  plot_data_reactome <- grouped_data_reactome %>%
    mutate(
      p.adjust = as.numeric(gsub(",", ".", p.adjust)),
      Count = as.numeric(Count)
    ) %>%
    filter(!is.na(p.adjust), !is.na(Theme))
  
  plot_data_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice(1:5) %>% # Show top 5 pathways per theme
    ungroup()
  
  theme_order_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    summarise(max_count = max(Count)) %>%
    arrange(desc(max_count)) %>%
    pull(Theme)
  
  plot_data_reactome <- plot_data_reactome %>%
    mutate(
      Description = reorder_within(Description, Count, Theme),
      Theme = factor(Theme, levels = rev(theme_order_reactome))
    )
  
  # --- Create the customized bar plot ---
  p_final_reactome <- ggplot(plot_data_reactome, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_y_reordered() +
    scale_fill_gradientn(colors = c("red", "magenta", "blue"), name = "p.adjust") +
    labs(title = "Unique Reactome Pathways for SIRT5_KO Downregulated Genes", x = "Count", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
      strip.background = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  # --- Display the final plot ---
  print(p_final_reactome)
}


# Extract the genes contributing to unique Reactome pathways 
# 1. Get the names of the pathways that are unique to SIRT5_KO
sirt5_unique_reactome_names <- unique_pathways_df_reactome %>%
  dplyr::filter(Cluster == "SIRT5_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df_reactome <- as.data.frame(compare_reactome_results)

# 3. Filter the full results
sirt5_unique_genes_df_reactome <- full_results_df_reactome %>%
  dplyr::filter(
    Cluster == "SIRT5_KO",
    Description %in% sirt5_unique_reactome_names
  )

# 4. Extract, split, and combine all the gene IDs
all_sirt5_unique_genes_reactome <- sirt5_unique_genes_df_reactome %>%
  pull(geneID) %>%
  strsplit(split = "/") %>%
  unlist() %>%
  unique()

# 5. Count them!
number_of_sirt5_genes_reactome <- length(all_sirt5_unique_genes_reactome)

cat("The number of genes contributing to the unique Reactome functions of SIRT5 is:",
    number_of_sirt5_genes_reactome, "\n")




######## For SIRT6
sirt6_unique_reactome_df <- list_of_unique_reactome_tables$SIRT6_KO

if (!is.null(sirt6_unique_reactome_df) && nrow(sirt6_unique_reactome_df) > 0) {
  
  # --- Export the table for manual thematic grouping ---
  write.csv(sirt6_unique_reactome_df, "SIRT6_unique_down_Reactome_for_grouping.csv", row.names = FALSE)
  
  cat("Table 'SIRT6_unique_down_Reactome_for_grouping.csv' was saved.\n")
  cat("Please open it in a spreadsheet editor, add a 'Theme' column, and save it as a new file.\n")
  
} else {
  cat("No unique Reactome pathways found for SIRT6_KO.\n")
}

# ! No unique Reactome pathways found for SIRT6_KO ! 




######## For SIRT7
sirt7_unique_reactome_df <- list_of_unique_reactome_tables$SIRT7_KO

if (!is.null(sirt7_unique_reactome_df) && nrow(sirt7_unique_reactome_df) > 0) {
  
  # --- Export the table for manual thematic grouping ---
  write.csv(sirt7_unique_reactome_df, "SIRT7_unique_down_Reactome_for_grouping.csv", row.names = FALSE)
  
  cat("Table 'SIRT7_unique_down_Reactome_for_grouping.csv' was saved.\n")
  cat("Please open it in a spreadsheet editor, add a 'Theme' column, and save it as a new file.\n")
  
} else {
  cat("No unique Reactome pathways found for SIRT7_KO.\n")
}

# Manually creat the "grouped" CSV file 

# Read your new, grouped table 
grouped_data_reactome <- read.csv2("SIRT7_unique_down_Reactome_grouped.csv", sep=";")

# Prepare the data for plotting 
if (nrow(grouped_data_reactome) > 0) {
  
  plot_data_reactome <- grouped_data_reactome %>%
    mutate(
      p.adjust = as.numeric(gsub(",", ".", p.adjust)),
      Count = as.numeric(Count)
    ) %>%
    filter(!is.na(p.adjust), !is.na(Theme))
  
  plot_data_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice(1:5) %>% # Show top 5 pathways per theme
    ungroup()
  
  theme_order_reactome <- plot_data_reactome %>%
    group_by(Theme) %>%
    summarise(max_count = max(Count)) %>%
    arrange(desc(max_count)) %>%
    pull(Theme)
  
  plot_data_reactome <- plot_data_reactome %>%
    mutate(
      Description = reorder_within(Description, Count, Theme),
      Theme = factor(Theme, levels = rev(theme_order_reactome))
    )
  
  # --- Create the customized bar plot ---
  p_final_reactome <- ggplot(plot_data_reactome, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    facet_grid(Theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_y_reordered() +
    scale_fill_gradientn(colors = c("red", "magenta", "blue"), name = "p.adjust") +
    labs(title = "Unique Reactome Pathways for SIRT7_KO Downregulated Genes", x = "Count", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
      strip.background = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  # --- Display the final plot ---
  print(p_final_reactome)
}


# Extract the genes contributing to unique Reactome pathways 
# 1. Get the names of the pathways that are unique to SIRT7_KO
sirt7_unique_reactome_names <- unique_pathways_df_reactome %>%
  dplyr::filter(Cluster == "SIRT7_KO") %>%
  pull(Description)

# 2. Go back to the original, full results data frame
full_results_df_reactome <- as.data.frame(compare_reactome_results)

# 3. Filter the full results
sirt7_unique_genes_df_reactome <- full_results_df_reactome %>%
  dplyr::filter(
    Cluster == "SIRT7_KO",
    Description %in% sirt7_unique_reactome_names
  )

# 4. Extract, split, and combine all the gene IDs
all_sirt7_unique_genes_reactome <- sirt7_unique_genes_df_reactome %>%
  pull(geneID) %>%
  strsplit(split = "/") %>%
  unlist() %>%
  unique()

# 5. Count them!
number_of_sirt7_genes_reactome <- length(all_sirt7_unique_genes_reactome)

cat("The number of genes contributing to the unique Reactome functions of SIRT7 is:",
    number_of_sirt7_genes_reactome, "\n")





### VIII. Upsetplot with common and unique up-/down-regulated DEGs ####

install.packages("UpSetR")
library(UpSetR)

# --- Create the input list for UPREGULATED genes ---
upregulated_gene_lists <- lapply(all_sig_tables, function(df) {
  # Filter for genes with log2FoldChange > 0 and extract the gene names
  df %>%
    dplyr::filter(log2FoldChange > 0) %>%
    pull(gene_name) %>%
    na.omit() # Remove any potential NA values
})


# --- Create the input list for DOWNREGULATED genes ---
downregulated_gene_lists <- lapply(all_sig_tables, function(df) {
  # Filter for genes with log2FoldChange < 0 and extract the gene names
  df %>%
    dplyr::filter(log2FoldChange < 0) %>%
    pull(gene_name) %>%
    na.omit()
})

# Let's check what we created
# It's a list where each element is a vector of gene names for that condition.
head(upregulated_gene_lists$SIRT1_KO)

# Define the desired order and colors for the sets on the left
set_order <- c("SIRT1_KO", "SIRT2_KO", "SIRT3_KO", "SIRT4_KO", "SIRT5_KO", "SIRT6_KO", "SIRT7_KO")
set_colors <- c("tan", "darkblue", "deepskyblue", "palevioletred4", "purple4", "darkslateblue", "darkcyan")


# CREATE A LIST OF QUERIES - ONE FOR EACH UNIQUE SET
# This is more robust than a single query for degree == 1.
# We will create one query for each sirtuin to highlight it with its specific color.
query_list <- lapply(seq_along(set_order), function(i) {
  list(
    # The query finds the intersection that corresponds ONLY to this set.
    query = intersects,
    params = list(set_order[i]), # The specific set, e.g., "SIRT1_KO"
    color = "deeppink4",
    active = TRUE
  )
})


# Create the Upset Plot for Upregulated Genes
upset(
  fromList(upregulated_gene_lists),
  
  # --- Core Parameters ---
  sets = rev(set_order),          # Enforce the set order
  nintersects = 24,              # Show the top 20 intersections
  order.by = "freq",             # Order by size
  
  # --- Key Customizations for Highlighting ---
  
  # Use the list of queries we just created.
  # The default color for bars not matching any query will be gray.
  queries = query_list,
  
  # Set the default color for bars that are NOT matched by a query (i.e., the shared sets)
  main.bar.color = "gray",
  
  # --- Other Aesthetic Parameters ---
  mainbar.y.label = "",
  sets.bar.color = rev(set_colors),
  text.scale = 1.5,
  set_size.show = FALSE
)


grid.text(
  "Unique upregulated DEGs (vs WT)",
  y = 0.95, # Располагаем вверху по центру (y=1 - самый верх)
  x = 0.72,
  gp = gpar(fontsize = 15, col = "deeppink4") # Устанавливаем размер и цвет
)

# Задача 1 (Часть 2): Добавляем метку оси Y ближе к оси
grid.text(
  "Gene number",
  x = 0.28, # Располагаем слева (x=0 - самый левый край)
  y = 0.7,
  rot = 90, # Поворачиваем на 90 градусов
  gp = gpar(fontsize = 13) # Устанавливаем размер шрифта
)





# Create the Upset Plot for Downregulated Genes
query_list <- lapply(seq_along(set_order), function(i) {
  list(
    # The query finds the intersection that corresponds ONLY to this set.
    query = intersects,
    params = list(set_order[i]), # The specific set, e.g., "SIRT1_KO"
    color = "blue",
    active = TRUE
  )
})


# Create the Upset Plot for Downregulated Genes
upset(
  fromList(downregulated_gene_lists),
  
  # --- Core Parameters ---
  sets = rev(set_order),          # Enforce the set order
  nintersects = 20,              # Show the top 20 intersections
  order.by = "freq",             # Order by size
  
  # --- Key Customizations for Highlighting ---
  
  # Use the list of queries we just created.
  # The default color for bars not matching any query will be gray.
  queries = query_list,
  
  # Set the default color for bars that are NOT matched by a query (i.e., the shared sets)
  main.bar.color = "gray",
  
  # --- Other Aesthetic Parameters ---
  mainbar.y.label = "",
  sets.bar.color = rev(set_colors),
  text.scale = 1.5,
  set_size.show = FALSE
)

grid.text(
  "Unique downregulated DEGs (vs WT)",
  y = 0.95, # Располагаем вверху по центру (y=1 - самый верх)
  x = 0.72,
  gp = gpar(fontsize = 15, col = "blue") # Устанавливаем размер и цвет
)

# Задача 1 (Часть 2): Добавляем метку оси Y ближе к оси
grid.text(
  "Gene number",
  x = 0.28, # Располагаем слева (x=0 - самый левый край)
  y = 0.7,
  rot = 90, # Поворачиваем на 90 градусов
  gp = gpar(fontsize = 13) # Устанавливаем размер шрифта
)






####################### 1. Extract unique up-/down-regulated DEGs #####

# Extract unique up/downregulated DEGs for the further enrichment analysis

# Now, loop through this list to find the unique genes for each condition
unique_upregulated_list <- list()
for (i in seq_along(upregulated_gene_lists)) {
  current_name <- names(upregulated_gene_lists)[i]
  current_genes <- upregulated_gene_lists[[i]]
  # Compare against the pool of all OTHER upregulated gene lists
  other_genes_pool <- unique(unlist(upregulated_gene_lists[-i]))
  
  unique_genes <- setdiff(current_genes, other_genes_pool)
  
  unique_upregulated_list[[current_name]] <- unique_genes
}


# Now, loop through this list to find the unique genes for each condition
unique_downregulated_list <- list()
for (i in seq_along(downregulated_gene_lists)) {
  current_name <- names(downregulated_gene_lists)[i]
  current_genes <- downregulated_gene_lists[[i]]
  # Compare against the pool of all OTHER downregulated gene lists
  other_genes_pool <- unique(unlist(downregulated_gene_lists[-i]))
  
  unique_genes <- setdiff(current_genes, other_genes_pool)
  
  unique_downregulated_list[[current_name]] <- unique_genes
}


# Check the number of unique upregulated genes for each condition
sapply(unique_upregulated_list, length)

# Check the number of unique downregulated genes for each condition
sapply(unique_downregulated_list, length)



######### 2. GO enrichment analysis on unique up-/down-regulated DEGs #####

# 1. For SIRT1_KO
sirt1_unique_up_genes <- unique_upregulated_list$SIRT1_KO
sirt1_unique_down_genes <- unique_downregulated_list$SIRT1_KO

# Run the enrichment analysis (with a safety check)
if (length(sirt1_unique_up_genes) > 5) { # Only run if you have a reasonable number of genes
  
  go_sirt1_unique_up <- enrichGO(
    gene          = sirt1_unique_up_genes,
    universe      = background, # Use your full background list
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 1,
    qvalueCutoff  = 0.5
  )
  
  # Visualize the results
  if (!is.null(go_sirt1_unique_up) && nrow(go_sirt1_unique_up) > 0) {
    dotplot(go_sirt1_unique_up, title = "GO for unique upregulated SIRT1_KO genes")
  } else {
    cat("No significant GO enrichment found for unique upregulated SIRT1_KO genes.\n")
  }
  
} else {
  cat("Not enough unique genes for SIRT1_KO to perform enrichment analysis.\n")
}

# Run the enrichment analysis (with a safety check)
if (length(sirt1_unique_down_genes) > 5) { # Only run if you have a reasonable number of genes
  
  go_sirt1_unique_down <- enrichGO(
    gene          = sirt1_unique_down_genes,
    universe      = background, # Use your full background list
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 1,
    qvalueCutoff  = 0.5
  )
  
  # Visualize the results
  if (!is.null(go_sirt1_unique_down) && nrow(go_sirt1_unique_down) > 0) {
    dotplot(go_sirt1_unique_down, title = "GO for unique downregulated SIRT1_KO genes")
  } else {
    cat("No significant GO enrichment found for unique downregulated SIRT1_KO genes.\n")
  }
  
} else {
  cat("Not enough unique genes for SIRT1_KO to perform enrichment analysis.\n")
}





# 2. For SIRT2_KO
sirt2_unique_up_genes <- unique_upregulated_list$SIRT2_KO
sirt2_unique_down_genes <- unique_downregulated_list$SIRT2_KO

# Run the enrichment analysis (with a safety check)
if (length(sirt2_unique_up_genes) > 5) { # Only run if you have a reasonable number of genes
  
  go_sirt2_unique_up <- enrichGO(
    gene          = sirt2_unique_up_genes,
    universe      = background, # Use your full background list
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  )
  
  # Visualize the results
  if (!is.null(go_sirt2_unique_up) && nrow(go_sirt2_unique_up) > 0) {
    dotplot(go_sirt2_unique_up, title = "GO for unique upregulated SIRT2_KO genes")
  } else {
    cat("No significant GO enrichment found for unique upregulated SIRT2_KO genes.\n")
  }
  
} else {
  cat("Not enough unique genes for SIRT2_KO to perform enrichment analysis.\n")
}

# Run the enrichment analysis (with a safety check)
if (length(sirt2_unique_down_genes) > 5) { # Only run if you have a reasonable number of genes
  
  go_sirt2_unique_down <- enrichGO(
    gene          = sirt2_unique_down_genes,
    universe      = background, # Use your full background list
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 1,
    qvalueCutoff  = 0.5
  )
  
  # Visualize the results
  if (!is.null(go_sirt2_unique_down) && nrow(go_sirt2_unique_down) > 0) {
    dotplot(go_sirt2_unique_down, title = "GO for unique downregulated SIRT2_KO genes")
  } else {
    cat("No significant GO enrichment found for unique downregulated SIRT2_KO genes.\n")
  }
  
} else {
  cat("Not enough unique genes for SIRT2_KO to perform enrichment analysis.\n")
}





# 3. For SIRT3_KO
sirt3_unique_up_genes <- unique_upregulated_list$SIRT3_KO
sirt3_unique_down_genes <- unique_downregulated_list$SIRT3_KO

# Run the enrichment analysis (with a safety check)
if (length(sirt3_unique_up_genes) > 5) { # Only run if you have a reasonable number of genes
  
  go_sirt3_unique_up <- enrichGO(
    gene          = sirt3_unique_up_genes,
    universe      = background, # Use your full background list
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 1,
    qvalueCutoff  = 0.5
  )
  
  # Visualize the results
  if (!is.null(go_sirt3_unique_up) && nrow(go_sirt3_unique_up) > 0) {
    dotplot(go_sirt3_unique_up, title = "GO for unique upregulated SIRT3_KO genes")
  } else {
    cat("No significant GO enrichment found for unique upregulated SIRT3_KO genes.\n")
  }
  
} else {
  cat("Not enough unique genes for SIRT3_KO to perform enrichment analysis.\n")
}

# Run the enrichment analysis (with a safety check)
if (length(sirt3_unique_down_genes) > 5) { # Only run if you have a reasonable number of genes
  
  go_sirt3_unique_down <- enrichGO(
    gene          = sirt3_unique_down_genes,
    universe      = background, # Use your full background list
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 1,
    qvalueCutoff  = 0.5
  )
  
  # Visualize the results
  if (!is.null(go_sirt3_unique_down) && nrow(go_sirt3_unique_down) > 0) {
    dotplot(go_sirt3_unique_down, title = "GO for unique downregulated SIRT3_KO genes")
  } else {
    cat("No significant GO enrichment found for unique downregulated SIRT3_KO genes.\n")
  }
  
} else {
  cat("Not enough unique genes for SIRT3_KO to perform enrichment analysis.\n")
}




# 4. For SIRT4_KO
sirt4_unique_up_genes <- unique_upregulated_list$SIRT4_KO
sirt4_unique_down_genes <- unique_downregulated_list$SIRT4_KO

# Run the enrichment analysis (with a safety check)
if (length(sirt4_unique_up_genes) > 5) { # Only run if you have a reasonable number of genes
  
  go_sirt4_unique_up <- enrichGO(
    gene          = sirt4_unique_up_genes,
    universe      = background, # Use your full background list
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 1,
    qvalueCutoff  = 0.5
  )
  
  # Visualize the results
  if (!is.null(go_sirt4_unique_up) && nrow(go_sirt4_unique_up) > 0) {
    dotplot(go_sirt4_unique_up, title = "GO for unique upregulated SIRT4_KO genes")
  } else {
    cat("No significant GO enrichment found for unique upregulated SIRT4_KO genes.\n")
  }
  
} else {
  cat("Not enough unique genes for SIRT4_KO to perform enrichment analysis.\n")
}

# Run the enrichment analysis (with a safety check)
if (length(sirt4_unique_down_genes) > 5) { # Only run if you have a reasonable number of genes
  
  go_sirt4_unique_down <- enrichGO(
    gene          = sirt4_unique_down_genes,
    universe      = background, # Use your full background list
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 1,
    qvalueCutoff  = 0.5
  )
  
  # Visualize the results
  if (!is.null(go_sirt4_unique_down) && nrow(go_sirt4_unique_down) > 0) {
    dotplot(go_sirt4_unique_down, title = "GO for unique downregulated SIRT4_KO genes")
  } else {
    cat("No significant GO enrichment found for unique downregulated SIRT4_KO genes.\n")
  }
  
} else {
  cat("Not enough unique genes for SIRT4_KO to perform enrichment analysis.\n")
}




# 5. For SIRT5_KO
sirt5_unique_up_genes <- unique_upregulated_list$SIRT5_KO
sirt5_unique_down_genes <- unique_downregulated_list$SIRT5_KO

# Run the enrichment analysis (with a safety check)
if (length(sirt5_unique_up_genes) > 5) { # Only run if you have a reasonable number of genes
  
  go_sirt5_unique_up <- enrichGO(
    gene          = sirt5_unique_up_genes,
    universe      = background, # Use your full background list
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 1,
    qvalueCutoff  = 0.5
  )
  
  # Visualize the results
  if (!is.null(go_sirt5_unique_up) && nrow(go_sirt5_unique_up) > 0) {
    dotplot(go_sirt5_unique_up, title = "GO for unique upregulated SIRT5_KO genes")
  } else {
    cat("No significant GO enrichment found for unique upregulated SIRT5_KO genes.\n")
  }
  
} else {
  cat("Not enough unique genes for SIRT5_KO to perform enrichment analysis.\n")
}

# Run the enrichment analysis (with a safety check)
if (length(sirt5_unique_down_genes) > 5) { # Only run if you have a reasonable number of genes
  
  go_sirt5_unique_down <- enrichGO(
    gene          = sirt5_unique_down_genes,
    universe      = background, # Use your full background list
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 1,
    qvalueCutoff  = 0.5
  )
  
  # Visualize the results
  if (!is.null(go_sirt5_unique_down) && nrow(go_sirt5_unique_down) > 0) {
    dotplot(go_sirt5_unique_down, title = "GO for unique downregulated SIRT5_KO genes")
  } else {
    cat("No significant GO enrichment found for unique downregulated SIRT5_KO genes.\n")
  }
  
} else {
  cat("Not enough unique genes for SIRT5_KO to perform enrichment analysis.\n")
}




# 6. For SIRT6_KO
sirt6_unique_up_genes <- unique_upregulated_list$SIRT6_KO
sirt6_unique_down_genes <- unique_downregulated_list$SIRT6_KO

# Run the enrichment analysis (with a safety check)
if (length(sirt6_unique_up_genes) > 5) { # Only run if you have a reasonable number of genes
  
  go_sirt6_unique_up <- enrichGO(
    gene          = sirt6_unique_up_genes,
    universe      = background, # Use your full background list
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 1,
    qvalueCutoff  = 0.5
  )
  
  # Visualize the results
  if (!is.null(go_sirt6_unique_up) && nrow(go_sirt6_unique_up) > 0) {
    dotplot(go_sirt6_unique_up, title = "GO for unique upregulated SIRT6_KO genes")
  } else {
    cat("No significant GO enrichment found for unique upregulated SIRT6_KO genes.\n")
  }
  
} else {
  cat("Not enough unique genes for SIRT6_KO to perform enrichment analysis.\n")
}

# Run the enrichment analysis (with a safety check)
if (length(sirt6_unique_down_genes) > 5) { # Only run if you have a reasonable number of genes
  
  go_sirt6_unique_down <- enrichGO(
    gene          = sirt6_unique_down_genes,
    universe      = background, # Use your full background list
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 1,
    qvalueCutoff  = 0.5
  )
  
  # Visualize the results
  if (!is.null(go_sirt6_unique_down) && nrow(go_sirt6_unique_down) > 0) {
    dotplot(go_sirt6_unique_down, title = "GO for unique downregulated SIRT6_KO genes")
  } else {
    cat("No significant GO enrichment found for unique downregulated SIRT6_KO genes.\n")
  }
  
} else {
  cat("Not enough unique genes for SIRT6_KO to perform enrichment analysis.\n")
}





# 7. For SIRT7_KO
sirt7_unique_up_genes <- unique_upregulated_list$SIRT7_KO
sirt7_unique_down_genes <- unique_downregulated_list$SIRT7_KO

# Run the enrichment analysis (with a safety check)
if (length(sirt7_unique_up_genes) > 5) { # Only run if you have a reasonable number of genes
  
  go_sirt7_unique_up <- enrichGO(
    gene          = sirt7_unique_up_genes,
    universe      = background, # Use your full background list
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 1,
    qvalueCutoff  = 0.5
  )
  
  # Visualize the results
  if (!is.null(go_sirt7_unique_up) && nrow(go_sirt7_unique_up) > 0) {
    dotplot(go_sirt7_unique_up, title = "GO for unique upregulated SIRT7_KO genes")
  } else {
    cat("No significant GO enrichment found for unique upregulated SIRT7_KO genes.\n")
  }
  
} else {
  cat("Not enough unique genes for SIRT7_KO to perform enrichment analysis.\n")
}

# Run the enrichment analysis (with a safety check)
if (length(sirt7_unique_down_genes) > 5) { # Only run if you have a reasonable number of genes
  
  go_sirt7_unique_down <- enrichGO(
    gene          = sirt7_unique_down_genes,
    universe      = background, # Use your full background list
    keyType       = "SYMBOL",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  )
  
  # Visualize the results
  if (!is.null(go_sirt7_unique_down) && nrow(go_sirt7_unique_down) > 0) {
    dotplot(go_sirt7_unique_down, title = "GO for unique downregulated SIRT7_KO genes")
  } else {
    cat("No significant GO enrichment found for unique downregulated SIRT7_KO genes.\n")
  }
  
} else {
  cat("Not enough unique genes for SIRT7_KO to perform enrichment analysis.\n")
}
