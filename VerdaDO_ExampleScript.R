library(DOSE)
library(clusterProfiler)
library(ggplot2)
library(tidyverse)
library(readr)
library(ChIPseeker)
library(regioneR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

setwd("/Users/scottt7/Desktop/")


# Load in all the data
Adrenal_clust_cs <- read_tsv("Adrenal_internalClusters_individualHMRs_cellspecific.txt", col_names = F)
colnames(Adrenal_clust_cs)<-c("Chr","Start","End","Length")
Adrenal_clust_cs <- Adrenal_clust_cs %>% 
  dplyr::select("Chr","Start","End") %>%
  mutate("ID"=row_number(), "Strand"=".", "OtherCol"="1")

Bcell_clust_cs <- read_tsv("Bcell_internalClusters_individualHMRs_cellspecific.txt", col_names = F)
colnames(Bcell_clust_cs)<-c("Chr","Start","End","Length")
Bcell_clust_cs <- Bcell_clust_cs %>% 
  dplyr::select("Chr","Start","End") %>%
  mutate("ID"=row_number(), "Strand"=".", "OtherCol"="1")


# Convert them to GRanges
Adrenal_clust_cs_granges <- toGRanges(as.data.frame(Adrenal_clust_cs))
Bcell_clust_cs_granges <- toGRanges(as.data.frame(Bcell_clust_cs))


Adrenal_clust_sh_granges <- toGRanges(as.data.frame(Adrenal_clust_sh))
Bcell_clust_sh_granges <- toGRanges(as.data.frame(Bcell_clust_sh))


# Combine into a list
list_of_regions <- list(Adrenal=Adrenal_clust_cs_granges,
                        Bcell=Bcell_clust_cs_granges,
                        Liver=Liver_clust_cs_granges,
                        fSpinal=fSpinal_clust_cs_granges,
                        H1ESC=H1ESC_clust_cs_granges,
                        Rvent=Rvent_clust_cs_granges)


# Apply AnnotatePeak to this list 
# GENE ASSIGNMENT
peakAnnoList <- lapply(list_of_regions, annotatePeak, TxDb=txdb,
                       tssRegion=c(-1000, 2000), verbose=FALSE)




# # Gene segment annotation
# # plotAnnoPie(peakAnno)
# # lapply(peakAnnoList, plotAnnoPie)
# plotAnnoBar(peakAnnoList)

# # Distance to TSs
# lapply(peakAnnoList, plotDistToTSS, title="Distance to Nearest TSS")
# plotDistToTSS(peakAnnoList, title="Distance to Nearest TSS")

############################


# Gene prep 
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)


# Venn Plot of genes
# vennplot(genes_5_noAdr)
# vennplot(genes_clust_sh_5_noAdr)

############### ENRICHMENT ANALYSES
# Don't know if you'll need this following line 
names(genes) = sub("_", "\n", names(genes))

compDO <- compareCluster(geneCluster   = genes,
                           fun           = "enrichDO", #enrichGO #enrichKEGG
                           pvalueCutoff  = 0.2,
                           pAdjustMethod = "BH")

dotplot(compDO, showCategory = 10, title = "DO Enrichment Analysis")

compDO_point1 <- compareCluster(geneCluster   = genes,
                           fun           = "enrichDO",
                           pvalueCutoff  = 0.1,
                           pAdjustMethod = "BH")
dotplot(compDO_point1, showCategory = 10, title = "DO Enrichment Analysis")


###### GO 
compGO <- compareCluster(geneCluster   = genes,
                              fun = "enrichGO",
                              OrgDb = org.Hs.eg.db,
                              ont = "BP",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.01,
                              qvalueCutoff  = 0.05)
dotplot(compGO, showCategory = 7, title = "GO Enrichment Analysis")


######### KEGG

compKEGG_point1 <- compareCluster(geneCluster   = genes,
                         fun           = "enrichKEGG",
                         pvalueCutoff  = 0.1,
                         pAdjustMethod = "BH")
dotplot(compKEGG_point1, showCategory = 10, title = "KEGG Enrichment Analysis")

compKEGG_point05 <- compareCluster(geneCluster   = genes,
                                   fun           = "enrichKEGG",
                                   pvalueCutoff  = 0.05,
                                   pAdjustMethod = "BH")
dotplot(compKEGG_point05, showCategory = 10, title = "KEGG Enrichment Analysis")

compKEGG_point01 <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.01,
                           pAdjustMethod = "BH")
dotplot(compKEGG_point01, showCategory = 10, title = "KEGG Enrichment Analysis")


# Wanted to know if the kidney overlap with Liver
# Was powered by liver-kidney genes identified by the Human Protein Atlas database
# Wanna look at the liver genes
as.data.frame(peakAnnoList$Liver)$geneId

Liver_gene_symbols <- unique(bitr(as.data.frame(peakAnnoList$Liver)$geneId, 
     fromType="ENTREZID",
     toType="SYMBOL",
     OrgDb = org.Hs.eg.db))


####### TESTs


# Try out ChipSeeker to annotate the nearest gene
# test_bed <- annotatePeak(Bcell_clust_cs, tssRegion=c(-1000,2000),
#                      TxDb=txdb, annoDb="org.Hs.eg.db")
# test <- annotatePeak(Bcell_clust_cs_granges, tssRegion=c(-1000,2000),
#              TxDb=txdb, annoDb="org.Hs.eg.db")
# 
# test_geneIDs <- as.data.frame(test)$geneId
# names(test_geneIDs) = sub("_", "\n", names(test_geneIDs))
# 
# testDO <- enrichDO(test_geneIDs,
#                    ont = "DO", 
#                    pvalueCutoff = 0.2, 
#                    pAdjustMethod = "BH",  
#                    minGSSize = 2, 
#                    qvalueCutoff = 0.2)
# dotplot(testDO, showCategory=15)
# 
# testGO <- enrichGO(test_geneIDs,
#                    ont = "BP", 
#                    pvalueCutoff = 0.2, 
#                    pAdjustMethod = "BH",  
#                    minGSSize = 2, 
#                    qvalueCutoff = 0.2,
#                    org.Hs.eg.db)
# 
# dotplot(testGO, showCategory=15)
