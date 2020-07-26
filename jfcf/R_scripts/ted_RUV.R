#
# https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html
# https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-3-rnaseq/rnaseq_gene_CSAMA2016.html
#

library(DESeq2)
library(RUVSeq)
library(tximport)
library(ensembldb)
library(tximportData)
library(EnsDb.Hsapiens.v75)

rhdf5::h5closeAll()
tx2gene <- transcripts(EnsDb.Hsapiens.v75, return.type="DataFrame")[, c("tx_id", "gene_id")]
colnames(tx2gene) <- c("TransID", "GeneID")

conts <- c("JFCF_6_T_1J_1_3C", "IIICF_a2", "VA13") # Controls
excel <- read.table("./samples/SAMPLES_Excel.tsv", sep="\t", header=TRUE)
excel <- excel[,c("Mortal", "Immortal")]
excel <- rbind(data.frame(Mortal="JFCF_6_T_1J_1_3C", Immortal="JFCF_6_T_1J_1_3C"), excel)
excel <- rbind(data.frame(Mortal="IIICF_a2", Immortal="IIICF_a2"), excel)
excel <- rbind(data.frame(Mortal="VA13", Immortal="VA13"), excel)

data <- read.table("./samples/SAMPLES_RNA.tsv", sep="\t", header=TRUE)
data <- data[with(data, order(Date)),]
data$DateSample <- as.factor(paste(data$Date, "_", data$Sample, sep=""))
data$File <- paste(paste("./RNA_data/", data$File, sep=""), "/abundance.h5", sep="")
data$Control <- FALSE
data[data$Sample %in% conts,]$Control <- TRUE
data$Date <- as.factor(data$Date)

stopifnot(all(file.exists(data$File)))
stopifnot(class(data$Sample) == "factor")

test <- function(row) {
  print(row)
  cont <- excel[row,] # What to test
  stopifnot(cont$Mortal %in% data$Sample)
  stopifnot(cont$Immortal %in% data$Sample)
  
  tmp <- data[data$Sample==as.character(cont$Mortal) | data$Sample==as.character(cont$Immortal),]
  
  has2017  <- "2017" %in% tmp$Date
  has2018  <- "2018" %in% tmp$Date
  hasBatch <- has2017 && has2018
  
  # Are we validating with controls?
  isTestControls <- all(tmp$Control)
  
  # Always RUV if we're validating or there're batch effects
  useRUV <- isTestControls || hasBatch
  
  tmp$Control <- F # This is our test, so no longer a control
  tmp <- rbind(tmp, data[data$Control & !(data$Sample %in% tmp$Sample),])
  
  #
  # We're going to make two tests, without and with RUV normalisation. For RUV, we'll use the controls as adjustment factors.
  #
  
  stopifnot(all(file.exists(tmp$File)))
  tmp1 <- tmp[!tmp$Control,] # No controls
  tmp2 <- tmp
  stopifnot(nrow(tmp1) == 6 || nrow(tmp1) == 9)
  
  txi1 <- tximport(tmp1$File, type="kallisto", tx2gene=tx2gene, ignoreTxVersion=TRUE) # No controls
  txi2 <- tximport(tmp2$File, type="kallisto", tx2gene=tx2gene, ignoreTxVersion=TRUE) # RUV controls
  colnames(txi1$counts) <- tmp1$ComputerName
  colnames(txi2$counts) <- tmp2$ComputerName
  
  tmp1$Sample <- factor(tmp1$Sample); tmp2$Sample <- factor(tmp2$Sample) # Make sure we drop unused levels
  tmp1$DateSample <- factor(tmp1$DateSample); tmp2$DateSample <- factor(tmp2$DateSample) # Make sure we drop unused levels
  tmp1$ComputerName <- factor(tmp1$ComputerName) # Make sure we drop the unused levels
  tmp2$ComputerName <- factor(tmp2$ComputerName) # Make sure we drop the unused levels
  stopifnot(class(tmp1$Date)=="factor" && class(tmp1$Sample)=="factor")
  stopifnot(class(tmp2$Date)=="factor" && class(tmp2$Sample)=="factor")
  
  mo  <- cont$Mortal # Mortal
  m17 <- paste("2017_", mo, sep="") # 2017 + Mortal
  m18 <- paste("2018_", mo, sep="") # 2018 + Mortal
  
  if (isTestControls) {
    tmp1$DateSample <- relevel(tmp1$DateSample, ref=as.character(m17)) # While validating, it's always 2017 vs 2018
    tmp2$DateSample <- relevel(tmp2$DateSample, ref=as.character(m17)) # While validating, it's always 2017 vs 2018
  } else {
    tmp1$Sample <- relevel(tmp1$Sample, ref=as.character(mo)) # Always test against mortal
    tmp2$Sample <- relevel(tmp2$Sample, ref=as.character(mo)) # Always test against mortal
  }
  
  prefix <- paste(cont$Mortal, "-", cont$Immortal, sep="") # Eg: VA13-VA13
  
  #
  # DESeq2 without RUV, differential analysis ignoring everything else
  #
  if (isTestControls) {
    dds1 <- DESeqDataSetFromTximport(txi1, tmp1, ~DateSample)
  } else {
    dds1 <- DESeqDataSetFromTximport(txi1, tmp1, ~Sample)
  }
  
  #
  # DESeq2 with controls by RUV. When we're testing controls, they all come from the same sample. Thus, we need to add "Date" to the samples. "DateSample"
  # is exactly what we need. Otherwise, we'll just test "Sample" but add latent factors (W_1) to the design formula.
  #
  
  x <- txi2$counts
  if (isTestControls) {
    # "DateSample" because they come from the same sample
    set.b <- newSeqExpressionSet(round(x), phenoData=data.frame(row.names=colnames(x), Date=tmp2$Date, DateSample=tmp2$DateSample))
  } else {
    set.b <- newSeqExpressionSet(round(x), phenoData=data.frame(row.names=colnames(x), Date=tmp2$Date, Sample=tmp2$Sample))
  }
  
  tmp2$Index <- c(1:nrow(tmp2)) # Make the next step easier
  
  cs <- unique(tmp2[tmp2$Control,]$Sample) # Controls
  cc <- c()
  
  for (i in c(1:length(cs))) {
    print(cs[i])
    cc <- c(cc, tmp2[tmp2$Sample == cs[i],]$Index)
  }
  
  scIdx <- matrix(cc, nrow=length(cs), byrow=T)
  
  # Apply RUVs normalisation
  set.a <- RUVSeq::RUVs(set.b, row.names(set.b), k=1, scIdx=scIdx)
  
  if (useRUV) {
    png(paste("11/", prefix, "_BeforePCA.png", sep="")); EDASeq::plotPCA(set.b); dev.off()
    png(paste("11/", prefix, "_AfterPCA.png", sep=""));  EDASeq::plotPCA(set.a); dev.off()
  }
  
  #
  # When testing for RUV, we always add the latent factor as the first variable in the formula. Multi-factor additive factor.
  #
  
  if (isTestControls) {
    dds2 <- DESeqDataSetFromMatrix(countData=counts(set.a), colData=pData(set.a), design=~W_1 + DateSample)
  } else {
    dds2 <- DESeqDataSetFromMatrix(countData=counts(set.a), colData=pData(set.a), design=~W_1 + Sample)
  }
  
  stopifnot(!is.null(dds1) && !is.null(dds2))
  dds1 <- DESeq(dds1); dds2 <- DESeq(dds2);
  
  c1 <- counts(dds1) # Counts for dds1
  c2 <- tmp1         # Meta-data for dds1
  
  m1 <- c1[,colnames(c1) %in% c2[c2$Sample %in% cont$Mortal,]$ComputerName]
  m2 <- c1[,colnames(c1) %in% c2[c2$Sample %in% cont$Immortal,]$ComputerName]
  if (!isTestControls) { stopifnot((ncol(m1) + ncol(m2)) == ncol(c1)) }
  
  a1 <- rowMeans(m1) # Mean for mortal samples (RUV correction has no absolute effect)
  a2 <- rowMeans(m2) # Mean for immortal samples (RUV correction has no absolute effect)
  s1 <- transform(m1, SD=apply(m1, 1, sd))$SD
  s2 <- transform(m2, SD=apply(m2, 1, sd))$SD
  
  if (isTestControls) {
    r1 <- results(dds1)
  } else {
    r1 <- results(dds1, contrast=c("Sample", as.character(cont$Immortal), as.character(cont$Mortal)))
  }
  
  if (isTestControls) {
    r2 <- results(dds2, contrast=c("DateSample", as.character(m17), as.character(m18))) # Does the controls work? (i.e. less differential genes?)
  } else {
    # Numerator level comes before denominator
    r2 <- results(dds2, contrast=c("Sample", as.character(cont$Immortal), as.character(cont$Mortal)))
  }
  
  stopifnot(all(row.names(r1) == names(a1)))
  stopifnot(all(row.names(r1) == names(a2)))    
  stopifnot(all(row.names(r2) == names(a1)))
  stopifnot(all(row.names(r2) == names(a2)))
  
  r1$MortalAverage <- a1; r1$ImmortallAverage <- a2; r1$MortalSD <- s1; r1$ImmortalSD <- s2
  r2$MortalAverage <- a1; r2$ImmortallAverage <- a2; r2$MortalSD <- s1; r2$ImmortalSD <- s2
  
  png(paste("11/", prefix, "_MA.png",  sep="")); DESeq2::plotMA(r1); dev.off()
  if (useRUV) { png(paste("11/", prefix, "_MA_RUV.png", sep="")); DESeq2::plotMA(r2); dev.off() }
  
  write.table(r1, file=paste("11/", prefix, ".tsv", sep=""),  sep="\t", quote=FALSE)
  if (useRUV) { write.table(r2, file=paste("11/", prefix, "_RUV.tsv", sep=""), sep="\t", quote=FALSE) }
  
  print(prefix)
}

for (i in 1:nrow(excel)) { test(i) }
#mclapply(1:nrow(excel), test, mc.cores=4)