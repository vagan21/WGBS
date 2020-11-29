# Import Libraries
library(DOSE)
library(clusterProfiler)

# Set working directory if you do that sort of thing
setwd("/Users/scottt7/Desktop/")

# Import data - Expects a list with separate gene names in rows 
simpmem.adrenal <- read.table("/Users/scottt7/Desktop/june19.all.4col.clusters.simplemembership.adrenal.annotatePeaks.genename.txt", header = FALSE, stringsAsFactors = default.stringsAsFactors())

# Transform data to vector
simpmem.adrenal.V <- as.vector(unlist(simpmem.adrenal))

# Covert gene symbols to both ENSEMBL and ENTREZ IDs
# Output list should have 4 columns (1 per ID convention)
simpmem.adrenal.IDmatrix <- bitr(simpmem.adrenal.V, fromType = "SYMBOL",
                            toType = c("ENSEMBL", "ENTREZID", "ENSEMBLPROT"),
                            OrgDb = "org.Hs.eg.db")






# BioMart
# Can use with ENSGXXXX.YY type
convert_geneID <- function(gene){
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  genes_1 <- getBM(filters="ensembl_gene_id_version",
    attributes= c("entrezgene_id"),
    values=gene,
    mart=mart)
  genes_2 <- drop_na(genes_1) %>% distinct()
  return(genes_2)
}

convert_geneID(LIST)