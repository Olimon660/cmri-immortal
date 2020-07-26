library(RUVSeq)
library(zebrafishRNASeq)
library(tximport)
library(ensembldb)
library(tximportData)
library(EnsDb.Hsapiens.v75)
library(preprocessCore)
library(RColorBrewer)
library(EDASeq)
library(DESeq2)
colors <- brewer.pal(3, "Set2")

listColumns(EnsDb.Hsapiens.v75)
tx2gene <- transcripts(EnsDb.Hsapiens.v75, return.type="DataFrame", columns = c("tx_id", "gene_name"))[, c("tx_id", "gene_name")]
colnames(tx2gene) <- c("TransID", "GeneName")

sample.excel = read.table("./samples/SAMPLES_Excel.tsv", sep = '\t', header = T, stringsAsFactors = F)
sample.rna = read.table("./samples/SAMPLES_RNA.tsv", sep = '\t', header = T, stringsAsFactors = F)


JF = sample.excel[grepl("JF", sample.excel$Immortal),c(2,3)]
ALT_JF = JF[JF$TMM == "ALT",]
TEL_JF = JF[JF$TMM == "TEL",]

sample.rna.JF.ALT = sample.rna[sample.rna$Sample %in% ALT_JF$Immortal,]
sample.rna.JF.TEL = sample.rna[sample.rna$Sample %in% TEL_JF$Immortal,]

#-----counts-------
count.mortal = c()
mortal.files = c("./RNA_data/JFCF6_1_CATBHANXX_CCGTCC/abundance.h5",
                 "./RNA_data/JFCF6_2_CATBHANXX_GTCCGC/abundance.h5",
                 "./RNA_data/JFCF6_3_CATBHANXX_GTGAAA/abundance.h5")
for (file in mortal.files) {
  tmp = tximport(file, type="kallisto", tx2gene=tx2gene, ignoreTxVersion=TRUE)
  count.mortal = cbind(count.mortal, tmp$counts)
}
colnames(count.mortal) = c("JFCF6_1_CATBHANXX_CCGTCC",
                           "JFCF6_2_CATBHANXX_GTCCGC",
                           "JFCF6_3_CATBHANXX_GTGAAA")
saveRDS(count.mortal, "./data/counts.genename.mortal.rds")

counts.alt = c()

for (i in 1:length(sample.rna.JF.ALT$File)){
  tmp = tximport(paste0("./RNA_data/",sample.rna.JF.ALT$File[i],"/abundance.h5"), 
                 type="kallisto", tx2gene=tx2gene, ignoreTxVersion=TRUE)
  counts.alt = cbind(counts.alt, tmp$counts)
}
colnames(counts.alt) = sample.rna.JF.ALT$File
saveRDS(counts.alt, "./data/counts.genename.alt.rds")

counts.tel = c()
for (i in 1:length(sample.rna.JF.TEL$File)){
  tmp = tximport(paste0("./RNA_data/",sample.rna.JF.TEL$File[i],"/abundance.h5"), 
                 type="kallisto", tx2gene=tx2gene, ignoreTxVersion=TRUE)
  counts.tel = cbind(counts.tel, tmp$counts)
}
colnames(counts.tel) = sample.rna.JF.TEL$File
saveRDS(counts.tel, "./data/counts.genename.tel.rds")

counts = cbind(counts.alt, counts.tel)
counts = as.data.frame(counts)

sample.names = c(sample.rna.JF.ALT$ComputerName, sample.rna.JF.TEL$ComputerName)
x <- as.factor(c(rep("ALT", 33), rep("TEL", 24)))
set <- newSeqExpressionSet(as.matrix(counts), phenoData = data.frame(x, row.names=colnames(counts)))
plotRLE(set, outline=FALSE, col=colors[x], las=2)
