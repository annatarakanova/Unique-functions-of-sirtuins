# Unique-functions-of-sirtuins
Characterization of the common and divergent transcriptional and functional signatures across the SIRT1-7 knockout models in human mesenchymal progenitor cells (hMSCs).

## The data
Bulk RNA-seq data from human mesenchymal progenitor cells (hMSCs) - WT vs SIRT1-7 KO (8 conditions) and 3 replicates per each condition.

The reference article: https://www.cell.com/developmental-cell/fulltext/S1534-5807(24)00107-2

## The goal
To characterize the unique functional signatures for each of the SIRT1-7 knockout model.

## The workflow
1. Preprocessing of bulk RNA-seq data (nf-core/rnaseq)
2. Differential expression analysis (DESeq2)
3. Split all DEGs into up- and down-regulated
4. Perform GO/Reactome enrichment analysis
5. Obtain only unique pathways
6. Group the pathways and find unique functions
7. Perform WGCNA (weighted gene correlation network analysis) to answer the question "Which group og genes behave similarly across all my samples, and do these groups relate to my conditions?"

## The results
Differential expression analysis (DESeq2) and functional interpretation showed differences in the transcription signatures of individual sirtuins, indicating their specific role in the regulation of cellular processes.

