## Enrichment analysis ---------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

setwd('/Users/osipova/Documents/LabDocs/Brood_parasites_analysis/')
# db = 'vidMac'
# db = 'vidCha'
db = 'indInd'
file_name = paste0('MK_test_', db, '_ncbi/gene.longest.mk.tsv')


## Load gene symbols and weights
file_df <- read.csv(file_name, header=TRUE, sep='\t')
gene_dos <- na.omit(file_df[c('gene', 'dos')])
gene_dos <- gene_dos[order(-gene_dos$dos), ]
gene_dos <- gene_dos[gene_dos$dos != 0, ]


## Perform GSEA with clusterProfiler package
gene_list <- gene_dos$dos
names(gene_list) = gene_dos$gene

gse <- gseGO(geneList         = gene_list,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = 'BH',
                pvalueCutoff  = 1,
                minGSSize = 35,
                maxGSSize = 500)

gse_result = gse@result
gsea_out_file = paste0('MK_test_', db, '_ncbi/gse.tsv')
write.table(gse_result[gse_result$pvalue < 0.01, ], row.names=F, quote=F, file=gsea_out_file, sep="\t")

# gseaplot(gse, geneSetID=14) ## vidMac


## Perform gene set enrichment analysis (clusterProfiler). DoS > 0
# pval_thresh = 0.05
pval_thresh = 0.1

for (dos in c('pos', 'neg')) {
  
  if (dos == 'pos') {
    genes = na.omit(file_df[(file_df$mk.raw.p.value < pval_thresh) & (file_df$dos > 0), ])$gene
  }
  else {
    genes = na.omit(file_df[(file_df$mk.raw.p.value < pval_thresh) & (file_df$dos < 0), ])$gene
  }
    
  enrich_res <- enrichGO(
    gene          = genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.5,
    pAdjustMethod = "BH",
    minGSSize     = 35,
    maxGSSize     = 500,
    universe = file_df$gene
  )
  
  enrich_result = enrich_res@result
  enrich_out_file = paste0('MK_test_', db, '_ncbi/', dos, '.enrichGO.tsv')
  write.table(enrich_result[enrich_result$pvalue < 0.01, ], row.names=F, quote=F, file=enrich_out_file, sep="\t")
}

 

