library("biomaRt")

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", path="/biomart/martservice")
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

data <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position"), mart=ensembl)
colnames(data) <- c("GeneID", "GeneName", "Chrom", "Start", "End")
data$GeneName <- make.unique(data$GeneName)

write.table(data, "data/EnsemblGRCh38Genes.csv", sep=",", row.names=F, quote=F)
