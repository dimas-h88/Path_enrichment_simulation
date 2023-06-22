library(limma)
library(GEOquery)
library(edgeR)
library(phantasus)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(data.table)
library(reactome.db)
library(fgsea)
library(aggregation)
library(msigdbr)
library(stats)
source('do_bi_fora.R')
source('fgsea_Multi.R')
source('camera_voom_m1.R')
source('camera_voom_m0.R')
source('fgsea_voom_m0.R')
source('fgsea_voom_m1.R')
source('fgseaResample.R')
## Creating summary tables to store data


summary_table_fgsea <- data.table(pathway = character(),
                                  pval = numeric(),
                                  padj = numeric(),
                                  seed = numeric(),
                                  N = numeric(),
                                  Method = character())


summary_table_camera <- data.table(pathway = character(),
                                   pval = numeric(),
                                   padj = numeric(),
                                   seed = numeric(),
                                   N = numeric(),
                                   Method = character())


summary_table_fora <- data.table (pathway = character(),
                                  pval = numeric(),
                                  padj = numeric(),
                                  seed = numeric(),
                                  N = numeric(),
                                  Method = character())

summary_table_SE <- data.table (pathway = character(),
                                  pval = numeric(),
                                  padj = numeric(),
                                  seed = numeric(),
                                  N = numeric(),
                                  Method = character())

# Loading data
raw_data <- read.gct("New_no_sym.gct")
# Taking summary count matrix
counts <-raw_data@assayData$exprs
entrez_ids <- select(org.Hs.eg.db, keys=rownames(counts), columns="ENTREZID", keytype="ENSEMBL",  uniqueRows=TRUE)
entrez_ids <- entrez_ids[!duplicated(entrez_ids$ENSEMBL),] #deleting multiple conversion
na_rows <- is.na(entrez_ids$ENTREZID) # deleting NA names
counts <- counts[!na_rows,]# deleting rows with NA ENTREZIDs
entrez_ids <- entrez_ids[!na_rows,]
rownames(counts) <- entrez_ids$ENTREZID # Changing names for further analysis
counts <- counts[complete.cases(rownames(counts)),]
counts <- counts[order(-rowSums(counts)),]
counts <- head(counts, 12000)

# Summary metadata
metadata <- raw_data@phenoData@data
colnames(metadata)[3] <- 'cell_type'

counts <- counts[,order(colnames(counts))]
metadata <- metadata[order(metadata$title),]


# Taking and filtering pathways by names
msigdbr_df <- msigdbr(species = 'Homo sapiens', category = 'C5', subcategory = 'GO:BP')
msigdbr_list = split(x = msigdbr_df$entrez_gene, f = msigdbr_df$gs_name)

# Filter the list to keep only the elements that have intersections with entrez_ids
pathways <- lapply(msigdbr_list, intersect, y = entrez_ids$ENTREZID)
pathways <- pathways[!is.na(pathways)]
pathways <- lapply(pathways, unique)
pathways <- lapply(pathways, as.character)


pathways <- pathways[lapply(pathways, length) >= 15]
pathways <- pathways[lapply(pathways, length) < 500]

for (n in c(4, 8, 10, 20, 40)) {
  for (i in 1:50) {
    
    print(n)
    print(i)
    sprintf('Starting M0 vs M0  with seen = %d and number of reps = %d', i, n)
    # Sampling
    set.seed(i)
    sample <- sample.int(48, n)
    samp_exprs <- counts[, sample]
    samp_meta <- metadata [sample, ]
    samp_meta$group <- rep(c("group1", "group2"), each = n/2)
    
    #differential expression analysis
    dds <- DESeqDataSetFromMatrix(countData = samp_exprs,
                                  colData =  samp_meta,
                                  design = ~ group)
    dds <- DESeq(dds)
    res <- results(dds)
    res <- res[!is.na(res$stat),]
    res <- res[!is.na(res$padj),]
    
    ## FGSEA
    result_table_fgsea <- fgsea_Multi(res, pathways)
    summary_table_fgsea <- rbind(summary_table_fgsea, result_table_fgsea)
    
    ## CAMERA
    result_table_camera <- camera_voom_m0(samp_exprs, samp_meta, dds)
    summary_table_camera <- rbind(summary_table_camera, result_table_camera)
    
    
    ## FORA
    foraRes <- do_bi_fora(res, pathways)
    
    result_table_fora <- data.table(pathway = foraRes$pathway,
                                    pval = foraRes$pval,
                                    padj = p.adjust(foraRes$pval, method = 'BH'),
                                    seed = i,
                                    N = n,
                                    Method = 'FORA')
    summary_table_fora <- rbind(summary_table_fora, result_table_fora)
    
    
    
    ## FGSEA_SE
    rwm <- res@rownames
    de <- data.table(as.data.frame(res), keep.rownames = TRUE)
    de <- de[, list(rn, 
                    log2FoldChange, lfcSE, 
                    pval=pvalue, stat)]
    
    stats <- de[, setNames(log2FoldChange, rn)]
    fgseaSERes <- fgseaMultilevelSE(pathways = pathways, 
                                    stats, se=de$lfcSE, 
                                    gseaParam=1, nResample = 11)
    
    result_table_SE <- data.table(pathway = fgseaSERes$pathway,
                                  pval = fgseaSERes$pval,
                                  padj = fgseaSERes$padj,
                                  seed = i,
                                  N = n,
                                  Method = 'FGSEA_SE')
    summary_table_SE <- rbind(summary_table_SE, result_table_SE)
    
    
  }
}
write.csv(summary_table_fora, './M0/summary_table_fora_m0.tsv') 
write.csv(summary_table_camera, './M0/summary_table_camera_m0.tsv') 
write.csv(summary_table_fgsea, './M0/summary_table_fgsea_m0.tsv')
write.csv(summary_table_SE, './M0/summary_table_fgsea_SE_m0.tsv')

