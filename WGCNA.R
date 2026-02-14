BiocManager::install(c("WGCNA", "impute"))
library(WGCNA)
library(DESeq2)
library(dplyr)
library(ggplot2)

allowWGCNAThreads() # Allow multi-threading for WGCNA, which will speed it up significantly.

###########################################################################################
######################     Prepare the input expression data     ##########################

# 1. Create DESeq object
# import raw feature count table 
data <- read.table("~/SIRT17_hMSCs/DESeq/salmon.merged.gene_counts.tsv", sep = "\t", header = TRUE)
head(data)

# filter of the low- and non-expressed genes (at least 10 counts per gene). Take only numeric columns 
data <- data[rowSums(data[, 3:ncol(data)]) > 10, ]

# Prepare the count matrix WITH gene IDs as rownames
count_matrix <- data[, 3:ncol(data)]
rownames(count_matrix) <- data$gene_id # Use ENSEMBL IDs as rownames

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

# Join two columns with _ and create a new column Factors
meta$Factors <- paste(meta$Genotype, meta$SampleID, sep = '_') 

# Create a new data type called factors. Factors are how R understands categorical data (i.e., groups).
# In statistical models like the one used by DESeq2, one of your groups must serve as the reference level or baseline. 
# All comparisons are made against this group.
meta$Factors <- factor(meta$Factors, levels = c("WT_WT", "KO_SIRT1", "KO_SIRT2", "KO_SIRT3", "KO_SIRT4", "KO_SIRT5", "KO_SIRT6", "KO_SIRT7"))

# Re-create the DESeqDataSet object
dds_factors <- DESeqDataSetFromMatrix(countData = round(count_matrix), 
                                      colData = meta, 
                                      design = ~ Factors)


# WGCNA requires a matrix of normalized and variance-stabilized expression data - rlog or vst tranformation
# 2. Perform a variance stabilizing transformation (vst is faster than rlog for many samples)
vsd <- vst(dds_factors, blind = TRUE)

# The WGCNA input needs to have samples as rows and genes as columns.
# 3. Transpose the assay() output.
data_transposed <- as.data.frame(t(assay(vsd)))

# 4. Check for genes and samples with too many missing values (not usually an issue with RNA-seq)
gsg <- goodSamplesGenes(data_transposed, verbose = 3)
gsg$allOK
# If the last statement returns TRUE, all genes have passed the cuts.
# If not, we run the following line of code to remove the offending genes and samples.
# if (!gsg$allOK) {
# datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
# }



###########################################################################################
######################     Choose the Soft-Thresholding Power (β)     #####################

# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(data_transposed, powerVector = powers, verbose = 5, networkType = "signed")

# Plot the results to find the best power
par(mfrow = c(1, 2))
cex1 <- 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
# Add a horizontal line at the recommended R^2 cutoff
abline(h = 0.85, col = "red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")

# --- DECISION ---
# Look at the first plot. You want to choose the LOWEST power where the line
# crosses the red line (R^2 > 0.85). Let's assume it's 6 for this example.
softPower <- 9



###########################################################################################
####################     Construct the Network and Identify Modules     ###################

# Build the network and find the gene modules in one go using blockwiseModules (constructs TOM - topological overlap matrix)
# This step can take a while...
net <- blockwiseModules(
  data_transposed,
  power = softPower,
  networkType = "signed", # Use "signed" for more specific modules
  TOMType = "signed",
  minModuleSize = 30,      # We want modules with at least 30 genes
  reassignThreshold = 0,
  mergeCutHeight = 0.25,   # A value of 0.25 merges modules that are 75% similar
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "SIRT_TOM",
  verbose = 3
)

# The 'net' object contains all the results. Let's see how many modules we found.
table(net$colors)
# Module 0 (grey) is for genes that didn't fit into any other module.
# Module is a cluster of genes that are highly co-expressed across your samples.

# Convert the numeric module labels to color names
moduleColors <- labels2colors(net$colors)

# Plot the Dendrogram with Module Colors
plotDendroAndColors(
  net$dendrograms[[1]],           # plot the dendrogram for the first block of genes (there're 4 blocks in total)
  moduleColors[net$blockGenes[[1]]], # The vector of colors for the genes in this block
  "Module colors",                # A label for the color bar
  dendroLabels = FALSE,           # We don't want to see individual gene labels
  hang = 0.03,                    # How far the dendrogram hangs down
  addGuide = TRUE,                # Adds a helpful guide line
  guideHang = 0.05,               # Position of the guide line
  main = "Cluster Dendrogram of Genes" # A main title for the plot
)




###########################################################################################
####################     Relate Modules to Experimental Traits     ########################

# Connect the unsupervised modules back to your SIRT knockouts
# --- 1. Prepare the Trait Data ---
# We need a data frame where rows are samples and columns are traits.
# Let's create binary traits: is the sample a SIRT1_KO? (1 for yes, 0 for no).
traitData <- as.data.frame(model.matrix(~ 0 + Factors, data = meta))

# Clean up the names
names(traitData) <- gsub("Factors", "", names(traitData))

# Check that sample names match
all.equal(rownames(data_transposed), rownames(traitData)) # Should be TRUE


# --- 2. Correlate Modules with Traits ---
# Get Module Eigengenes (the summary profile for each module/hypothetical central gene))
MEs0 <- moduleEigengenes(data_transposed, net$colors)$eigengenes
# Reorder modules so similar modules are next to each other
MEs <- orderMEs(MEs0)

# Calculate the correlation and p-values
moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(data_transposed))

# --- Step 3: Convert the Correlation and P-value Matrices into a "Tidy" Data Frame ---
# This is the key step to prepare the data for ggplot2.

# Convert the correlation matrix to a data frame
cor_df <- as.data.frame(moduleTraitCor)
cor_df$module <- rownames(cor_df)

# Convert the p-value matrix to a data frame
pval_df <- as.data.frame(moduleTraitPvalue)
pval_df$module <- rownames(pval_df)

library(tidyr) # For the pivot_longer function

# "Melt" or "pivot" the data frames from wide to long format
cor_long <- cor_df %>%
  pivot_longer(
    cols = -module,
    names_to = "trait",
    values_to = "correlation"
  )

pval_long <- pval_df %>%
  pivot_longer(
    cols = -module,
    names_to = "trait",
    values_to = "p_value"
  )

# Join the two long data frames into one master data frame for plotting
plot_df <- left_join(cor_long, pval_long, by = c("module", "trait"))


# --- Step 4: Create the Heatmap using ggplot2 ---
# Create the plot
p_heatmap <- ggplot(plot_df, aes(x = trait, y = module, fill = correlation)) +
  # Use geom_tile to create the heatmap squares
  geom_tile(color = "white") +
  
  # --- ADD THIS MODIFIED geom_text() LAYER ---
  geom_text(
    aes(label = ifelse(
      # --- This is the new, complex condition ---
      abs(correlation) >= 0.70,
      
      # --- The rest is the same ---
      round(correlation, 2), # Value if TRUE
      ""                     # Value if FALSE
    )
    ),
    color = "black",
    size = 3
  ) +
  
  # --- The rest of the code is the same ---
  
  # Define the color scale
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1, 1),
    name = "Correlation"
  ) +
  
  # Apply a clean theme and customize
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  
  # Add a title
  labs(title = "Module-Trait Relationship") +
  
  # Ensure the aspect ratio is correct
  coord_fixed()

# --- Step 5: Display the Final Plot ---
print(p_heatmap)

# What does the heatmap show? - "Which of my co-expressed gene 'cliques' (modules) are associated with which of my experimental conditions (SIRT knockouts)?"



###########################################################################################
##############################     Enrichment analysis     ################################

# Let's perform enrichment analysis for module genes that I filtered (visualize on the heatmap)

# Step 1. Define modules of interest:
# Upregulated modules
sirt1_up_modules <- c("ME10")
sirt2_up_modules <- c("ME21", "ME27", "ME28")
sirt6_up_modules <- c("ME7")
sirt7_up_modules <- c("ME23")
wt_up_modules    <- c("ME1")

# Downregulated modules
sirt1_down_modules <- c("ME9")
sirt2_down_modules <- c("ME19")
sirt3_down_modules <- c("ME24")
sirt6_down_modules <- c("ME3")
sirt7_down_modules <- c("ME13")
wt_down_modules    <- c("ME2", "ME6", "ME17")


# Step 2. Upload libraries 
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)



# Step 3. Create a Map from ME Name to Module Color
# Get the numeric module labels for each gene from the 'net' object
module_numbers <- net$colors # This is a vector of numbers: 0, 1, 2, 3...
# Convert these numbers to their corresponding color names
module_colors <- labels2colors(module_numbers)

# Create the map: ME1 -> color, ME2 -> color, etc.
# The names of the MEs object are ME1, ME2, ME24, etc.
me_to_color_map <- data.frame(
  me_name = names(MEs), # e.g., "ME24", "ME4", "ME19"
  # The color is derived from the numeric part of the ME name
  color = labels2colors(as.numeric(substring(names(MEs), 3)))
)
print("--- ME Name to Color Map ---")
print(head(me_to_color_map))




# Step 4. Make a function that extract genes and perform enrichment analysis
perform_module_GO_analysis <- function(me_names_from_table, me_map, all_genes_and_colors, universe, analysis_title) {
  
  cat(paste("\n--- Running GO Analysis for:", analysis_title, "---\n"))
  
  # Find the corresponding colors for our ME names
  target_colors <- me_map$color[me_map$me_name %in% me_names_from_table]
  
  if (length(target_colors) == 0) {
    cat("No matching colors found for the given ME names. Skipping.\n")
    return(NULL)
  }
  
  # Extract all genes belonging to these colors.
  # These are your Gene Symbols.
  module_gene_symbols <- names(all_genes_and_colors)[all_genes_and_colors %in% target_colors]
  
  cat("Number of Gene Symbols extracted from module(s):", length(module_gene_symbols), "\n")
  
  # --- THIS IS THE CRUCIAL CONVERSION STEP (SYMBOL -> ENTREZID) ---
  
  # Use bitr() to translate GENE SYMBOL to ENTREZID
  # Note the change in fromType
  id_conversion_df <- bitr(module_gene_symbols,
                           fromType = "SYMBOL",      # <-- The type of ID we ARE GIVING
                           toType = "ENTREZID",    # The type of ID we WANT
                           OrgDb = org.Hs.eg.db)
  
  # Extract the clean vector of Entrez IDs
  module_entrez_ids <- id_conversion_df$ENTREZID
  
  cat("Number of genes successfully mapped to Entrez IDs:", length(module_entrez_ids), "\n")
  
  # --- End of conversion step ---
  
  # Perform enrichment analysis using the NEW Entrez ID list
  if (length(module_entrez_ids) > 10) {
    
    go_result <- enrichGO(
      gene = module_entrez_ids,   # <-- Use the converted IDs
      universe = universe,         # Your background MUST also be Entrez IDs
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",      # This now correctly matches our input
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.1
    )
    
    # (The rest of the plotting code is the same...)
    if (!is.null(go_result) && nrow(go_result) > 0) {
      go_result_simple <- simplify(go_result, cutoff = 0.7, by = "p.adjust")
      p <- dotplot(go_result_simple, showCategory = 10) + ggtitle(analysis_title)
      print(p)
      return(go_result)
    } else {
      cat("No GO enrichment found.\n")
      return(NULL)
    }
  } else {
    cat("Not enough genes in the module to perform GO analysis.\n")
    return(NULL)
  }
}



# Step 5. Run the GO analysis for each group from your table
# Create the vector of all gene names and their assigned colors
# (Assuming 'module_colors' and 'datExpr0' are from your WGCNA run)
all_genes_and_colors <- labels2colors(net$colors)
names(all_genes_and_colors) <- colnames(data_transposed)

# The background/universe should be a vector of all Entrez IDs tested by DESeq2
# (Assuming you have this prepared as 'background_entrez_ids')

# Run for SIRT1-specific upregulated modules
sirt1_up_go <- perform_module_GO_analysis(
  me_names_from_table = sirt1_up_modules, # e.g., c("ME10")
  me_map = me_to_color_map,
  all_genes_and_colors = all_genes_and_colors,
  universe = background_entrez_ids, # from DEG R script 
  analysis_title = "GO for SIRT1 upregulated modules (ME10)"
)

# Run for SIRT1-specific downregulated modules
sirt1_down_go <- perform_module_GO_analysis(
  me_names_from_table = sirt1_down_modules, 
  me_map = me_to_color_map,
  all_genes_and_colors = all_genes_and_colors,
  universe = background_entrez_ids, # from DEG R script 
  analysis_title = "GO for SIRT1 downregulated modules (ME9)"
)

# Run for SIRT2-specific upregulated modules
sirt2_up_go <- perform_module_GO_analysis(
  me_names_from_table = sirt2_up_modules, 
  me_map = me_to_color_map,
  all_genes_and_colors = all_genes_and_colors,
  universe = background_entrez_ids, # from DEG R script 
  analysis_title = "GO for SIRT2 upregulated modules (ME21, ME27, ME28)"
)

# Run for SIRT2-specific downregulated modules
sirt2_down_go <- perform_module_GO_analysis(
  me_names_from_table = sirt2_down_modules, 
  me_map = me_to_color_map,
  all_genes_and_colors = all_genes_and_colors,
  universe = background_entrez_ids, # from DEG R script 
  analysis_title = "GO for SIRT2 downregulated modules (ME19)"
)

# Run for SIRT3-specific downregulated modules
sirt3_down_go <- perform_module_GO_analysis(
  me_names_from_table = sirt3_down_modules, 
  me_map = me_to_color_map,
  all_genes_and_colors = all_genes_and_colors,
  universe = background_entrez_ids, # from DEG R script 
  analysis_title = "GO for SIRT3 downregulated modules (ME24)"
)

# Run for SIRT6-specific upregulated modules
sirt6_up_go <- perform_module_GO_analysis(
  me_names_from_table = sirt6_up_modules, 
  me_map = me_to_color_map,
  all_genes_and_colors = all_genes_and_colors,
  universe = background_entrez_ids, # from DEG R script 
  analysis_title = "GO for SIRT6 upregulated modules (ME7)"
)

# Run for SIRT6-specific downregulated modules
sirt6_down_go <- perform_module_GO_analysis(
  me_names_from_table = sirt6_down_modules, 
  me_map = me_to_color_map,
  all_genes_and_colors = all_genes_and_colors,
  universe = background_entrez_ids, # from DEG R script 
  analysis_title = "GO for SIRT6 downregulated modules (ME3)"
)

# Run for SIRT7-specific upregulated modules
sirt7_up_go <- perform_module_GO_analysis(
  me_names_from_table = sirt7_up_modules, 
  me_map = me_to_color_map,
  all_genes_and_colors = all_genes_and_colors,
  universe = background_entrez_ids, # from DEG R script 
  analysis_title = "GO for SIRT7 upregulated modules (ME23)"
)

# Run for SIRT7-specific downregulated modules
sirt7_down_go <- perform_module_GO_analysis(
  me_names_from_table = sirt7_down_modules, 
  me_map = me_to_color_map,
  all_genes_and_colors = all_genes_and_colors,
  universe = background_entrez_ids, # from DEG R script 
  analysis_title = "GO for SIRT7 downregulated modules (ME13)"
)

# Run for WT-specific upregulated modules
wt_up_go <- perform_module_GO_analysis(
  me_names_from_table = wt_up_modules,
  me_map = me_to_color_map,
  all_genes_and_colors = all_genes_and_colors,
  universe = background_entrez_ids,
  analysis_title = "GO for WT upregulated modules (ME1)"
)

# Run for WT-specific downregulated modules
wt_down_go <- perform_module_GO_analysis(
  me_names_from_table = wt_down_modules,
  me_map = me_to_color_map,
  all_genes_and_colors = all_genes_and_colors,
  universe = background_entrez_ids,
  analysis_title = "GO for WT downregulated modules (ME2, ME6, ME17)"
)
