##Project: Access Array Soay Sheep
##Aim: Sequencing pre-processing and denoising, merge by primer combination.
n##Author: Emanuel Heitlinger mod. by Víctor Hugo Jarquín-Díaz

##Call required libraries library(MultiAmplicon,
## lib.loc="/usr/local/lib/R/site-library/") using the local
## (devolopment version, if that branch is checked out) of
## MultiAmplicon
devtools::load_all("../MultiAmplicon")

library(ggplot2)
library(dada2)
library(reshape)
library(phyloseq)
library(data.table)
library(taxonomizr)
library(taxize)
library(parallel)
library(pheatmap)
library(tidyr)
library(dplyr)
library(ShortRead)

## re-run or use pre-computed results for different parts of the pipeline:
## Set to FALSE to use pre-computed and saved results, TRUE to redo analyses.
doQualEval <- FALSE

doFilter <- FALSE

doMultiAmpSort <- TRUE

doMultiAmpError <- TRUE

doMultiAmpPipe <- TRUE

doTax <- TRUE

doPS<- TRUE

###################dada2 pipeline#######################

path <- c(
  ## Soay Sheep (Multiamplicon run) Full sequencing Run (Good run)
  "2021_16_soay_Main_Run", 
  ## Hyena Pool 1 (Multiamplicon run) Preliminary test sequencing Run
  "2021_16_soay_Test_Run") 


fullpath <- paste0("/SAN/Victors_playground/AA_Soay/sequences/", path)

names(fullpath) <- path

fastqList <- lapply(fullpath, function (path) { 
  fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) 
  fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE)
  fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE)
  list(fastqF=fastqF, fastqR=fastqR)
})

if(doQualEval){
  readlenght <- lapply(fastqList, function (x) {
    con <- file(x[["fastqF"]][1],"r")
    ## first line
    secondLine <- readLines(con, n=2)[[2]]
    ### simple check whether it's V2 or V3 data
    nchar(secondLine)
  })
  
  allFastqF <- lapply(fastqList, function (x) {
    readFastq(x[["fastqF"]])
  })
  
  allFastqR <- lapply(fastqList, function (x) {
    readFastq(x[["fastqR"]])
  })
  
  
  sampleQual <- function (x) {
    ## sample quality scores of 100,000 sequences 
    qmat <- as(quality(x)[sample(100000)], "matrix")
    cols <- seq(1, ncol(qmat), by=10)
    sapply(cols, function (i) {
      mean(qmat[, i], na.rm=TRUE)
    })
  }
  
  qualityF <- lapply(allFastqF, sampleQual)
  qualityR <- lapply(allFastqR, sampleQual)
  
  shouldL <- max(unlist(lapply(qualityF, length)))
  
  qualityFilledF <- lapply(qualityF, function (x) {
    c(x, rep(NA, times=shouldL - length(x)))
  })
  
  qualityFilledR <- lapply(qualityR, function (x) {
    c(x, rep(NA, times=(shouldL - length(x))))
  })
  
  
  qualityDFF <- Reduce("cbind",  qualityFilledF)
  qualityDFR <- Reduce("cbind",  qualityFilledR)
  
  colnames(qualityDFF) <- path
  colnames(qualityDFR) <- path
  
  qualityDFFL <- reshape2::melt(qualityDFF)
  qualityDFFL$direction <- "forward"
  
  qualityDFRL <- reshape2::melt(qualityDFR)
  qualityDFRL$direction <- "reverse"
  
  qualityDFL <- rbind(qualityDFFL, qualityDFRL)
  
  qualityDFL$position <- qualityDFL$Var1*10 -10
  
  ggplot(qualityDFL, aes(position, value, color=Var2)) +
    geom_line() +
    facet_wrap(~direction)
}

## concluding from this that we can truncate: 
# at 200 for Rev # at 240 for Fwd

samplesList <- lapply (fastqList, function (x){
  samples<- gsub("-", "_", basename(x[["fastqF"]]))
  samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", samples)
  samples<- gsub("S\\d+_", "\\1", samples)
  paste(basename(dirname(x[["fastqF"]])), samples, sep="_")
})

fastqFall <- unlist(lapply(fastqList, "[[", "fastqF"))
fastqRall <- unlist(lapply(fastqList, "[[", "fastqR"))

samplesAll <- unlist(samplesList)

#Creation of a folder for filtrated reads 
filt_path <- "/SAN/Victors_playground/AA_Soay/tmp/filtered_all"

if(!file_test("-d", filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(samplesAll, "_F_filt.fastq.gz"))
names(filtFs) <- samplesAll
filtRs <- file.path(filt_path, paste0(samplesAll, "_R_filt.fastq.gz"))
names(filtRs) <- samplesAll

if(doFilter){
  filter.track <- lapply(seq_along(fastqFall),  function (i) {
    filterAndTrim(fastqFall[i], filtFs[i], fastqRall[i], filtRs[i],
                  truncLen=c(240,200), minLen=c(240,200), 
                  maxN=0, maxEE=2, truncQ=2, 
                  compress=TRUE, verbose=TRUE, multithread=TRUE,
                  matchIDs=TRUE) ## forward and reverse not matching otherwise 
  })
  saveRDS(filter.track, file="/SAN/Victors_playground/AA_Soay/tmp/filter.Rds")
} else {
  filter.track <- readRDS(file="/SAN/Victors_playground/AA_Soay/tmp/filter.Rds")
}

##Check the proportion of reads that passed the filtering 
filter <- as.data.frame(do.call(rbind, filter.track))
sum(filter[,"reads.out"])/sum(filter[,"reads.in"])

### Over 70% passed for all runs...
filter$run <- unlist(lapply(strsplit(samplesAll, "_Run_"), "[", 1))
##filter$run <- gsub("\\d_part\\d", "", filter$run)

## 76% in the main run and 61% in the test run
by(filter, filter$run, function (x) sum(x[,"reads.out"]/sum(x[,"reads.in"])))


files <- PairedReadFileSet(filtFs, filtRs)

### Get indexing data and indicate replicates
sampleIDs <- read.csv("/SAN/Victors_playground/AA_Soay/data/Soay_Sheep_Sample_Indexing.csv")
sampleIDs%>%
  select(Sequencing_ID)-> sampleIDs
sampleIDs <- as.data.frame(sampleIDs)

colnames(sampleIDs)[colnames(sampleIDs)%in%"Sequencing_ID"] <- "sampleID"

filter$sampleID <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", rownames(filter))
filter$sampleID<- gsub("-", "_", filter$sampleID)
filter$sampleID <- gsub("S\\d+_", "\\1", filter$sampleID)

filter$SnumIDs <- gsub("(S\\d{3,4})\\.(P\\d)\\.(FLD\\d{4}).*", "\\1_\\2_\\3",
                       rownames(filter))

sampleIDs <- merge(sampleIDs, filter, by="sampleID", all=TRUE)

sampleIDs$replicate <- ifelse(grepl("P1", sampleIDs$sampleID),
                              "Replicate_1", ifelse(grepl("P2", sampleIDs$sampleID), "Replicate_2", NA))

###################MultiAmplicon pipeline#######################
#Preparation of primer file
#Primers used in the arrays
ptable <- read.csv(file = "Data/primer.file.csv",
                   sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "Seq_F"]
primerR <- ptable[, "Seq_R"]
names(primerF) <- as.character(ptable[, "Name_F"])
names(primerR) <- as.character(ptable[, "Name_R"])
primer <- PrimerPairsSet(primerF, primerR)

M1 <- MultiAmplicon(primer, files)

sampleIDs$run<- paste(sampleIDs$run, "Run", sep="_")
rownames(sampleIDs) <- make.unique(paste(sampleIDs$run,
                                         gsub("-", "_", sampleIDs$sampleID),
                                         sep="_"))
##Control that both rownames match 
setdiff(rownames(sampleIDs), rownames(M1@sampleData))

MA <- addSampleData(M1, sampleIDs)

##Multi amplicon pipeline
if(doMultiAmpSort){
  filedir <- "/SAN/Victors_playground/AA_Soay/tmp/stratified_all"
  if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
  ## This step sort the reads into amplicons based on the number of primer pairs
  MA <- sortAmplicons(MA, n=1e+07, filedir=filedir) 
  ## pdf("Figures/overview_all_heat.pdf", width=16, height=31)
  ## pheatmap(log10(getRawCounts(MA)+1), show_colnames = FALSE)
  ## dev.off()
  saveRDS(MA, file="/SAN/Victors_playground/AA_Soay/tmp/MA_sorted.Rds")
} else {
  MA <- readRDS(file="/SAN/Victors_playground/AA_Soay/tmp/MA_sorted.Rds")
}

## separate sample data for each run

if(doMultiAmpError){
    MAF <- MA[, getSampleData(MA)$run%in%"2021_16_soay_Main_Run"]
    ## doing this only for the final run for now
    MAR <- derepMulti(MAF, mc.cores = 12)
    MAD <- dadaMulti(MAR, Ferr=NULL, selfConsist=TRUE,
                     Rerr=NULL,  pool=TRUE,
                     verbose=0, mc.cores = 12)
    saveRDS(MAD, file="/SAN/Victors_playground/AA_Soay/tmp/MAD.Rds")
} else {
    MAD <- readRDS(file="/SAN/Victors_playground/AA_Soay/tmp/MAD.Rds")
}

if(doMultiAmpPipe){
  MAM <- mergeMulti(MAD, mc.cores=12)
  propMerged <- calcPropMerged(MAM)
  
  MAM <- mergeMulti(MAM, justConcatenate=propMerged<0.7)
  MAS <- makeSequenceTableMulti(MAM, mc.cores=12)
  MA <- removeChimeraMulti(MAS, mc.cores=12) 
  
  ## Removing previous MA that consume a lot of space in RAM
  ##  rm(MAM, MAS, MAD, MAR, MAF)
  saveRDS(MA, file="/SAN/Victors_playground/AA_Soay/tmp/MA_piped.Rds")
} else {
  MA <- readRDS(file="/SAN/Victors_playground/AA_Soay/tmp/MA_piped.Rds")
}

rownames(MA@sampleData) <- getSampleData(MA)$sampleID

## AData from prepare_samples.R <- SOURCE ME!!!
M <- addSampleData(MA, AData)

## trackingF <- getPipelineSummary(MA) 
## PipSum <- plotPipelineSummary(trackingF) + scale_y_log10()
## ggsave("/SAN/Victors_playground/AA_Soay/tmp/Pipeline_track_1.pdf", PipSum,height = 15, width = 15) ##Temporal storage



##Lets run BLAST
if (doTax) {
  MA.1 <- blastTaxAnnot(MA.1,
                        db = "/SAN/db/blastdb/nt/nt",
                        negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                        infasta = "/SAN/Victors_playground/AA_Soay/output/Soay_1_in.fasta",
                        outblast = "/SAN/Victors_playground/AA_Soay/output/Soay_1_out.blt",
                        taxonSQL = "/SAN/db/taxonomy/taxonomizr.sql", 
                        num_threads = 20)

  MA.2 <- blastTaxAnnot(MA.2,
                        db = "/SAN/db/blastdb/nt/nt",
                        negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                        infasta = "/SAN/Victors_playground/AA_Soay/output/Soay_2_in.fasta",
                        outblast = "/SAN/Victors_playground/AA_Soay/output/Soay_2_out.blt",
                        taxonSQL = "/SAN/db/taxonomy/taxonomizr.sql", 
                        num_threads = 20)  
  saveRDS(MA.1, file="/SAN/Victors_playground/AA_Soay/output/MultiAmpliconTaxa_Soay_1.Rds") ###Results from full run 
  saveRDS(MA.2, file="/SAN/Victors_playground/AA_Soay/output/MultiAmpliconTaxa_Soay_2.Rds") ###Results from test run 
  
  plotAmpliconNumbers(MA.1)
  plotAmpliconNumbers(MA.2)
}else{
  MA.1<- readRDS(file="/SAN/Victors_playground/AA_Soay/output/MultiAmpliconTaxa_Soay_1.Rds")
  MA.2<- readRDS(file="/SAN/Victors_playground/AA_Soay/output/MultiAmpliconTaxa_Soay_2.Rds")   
  }

if (doPS) {
  ##Create individual Phyloseq objects 
  PS.1 <- toPhyloseq(MA.1, samples=colnames(MA.1))
  PS.2 <- toPhyloseq(MA.2, samples=colnames(MA.2))
  
  ##Separated list of primers
  PS.1.l <- toPhyloseq(MA.1, samples=colnames(MA.1), multi2Single = F) 
  
  ##PS.2.l <- toPhyloseq(MA.2, samples=colnames(MA.2), multi2Single = F)
  ##Error in access(object, "otu_table", errorIfNULL) : 
  ##otu_table slot is empty.
  
  ##Merge Phyloseq objects into a single 
  PS <- merge_phyloseq(PS.1, PS.2) 
  
  ##Store it for further analysis
  saveRDS(PS, file="/SAN/Victors_playground/AA_Soay/output/PhyloSeqCombi_Soay.Rds") ###Results from full + test run 
}

###Get rowcounts by primer pair 
rawcounts.1 <- rowSums(getRawCounts(MA.1))
rawcounts.1 <- data.frame(rawcounts.1)
rawcounts.1[,2] <- rownames(rawcounts.1)
colnames(rawcounts.1) <- c("Raw_counts_Main", "Primer_name")
rownames(rawcounts.1) <- c(1:nrow(rawcounts.1))
rawcounts.1 <- data.frame(Primer_name = rawcounts.1$Primer_name, Raw_counts_Main = rawcounts.1$Raw_counts_Main) ###change the order of the columns
rawcounts.1$Primer_name <- gsub(pattern = " ", replacement = "", x = rawcounts.1$Primer_name)
rawcounts.1$Primer_name <- gsub(pattern = "-", replacement = "_", x = rawcounts.1$Primer_name)

rawcounts.2 <- rowSums(getRawCounts(MA.2))
rawcounts.2 <- data.frame(rawcounts.2)
rawcounts.2[,2] <- rownames(rawcounts.2)
colnames(rawcounts.2) <- c("Raw_counts_Test", "Primer_name")
rownames(rawcounts.2) <- c(1:nrow(rawcounts.2))
rawcounts.2 <- data.frame(Primer_name = rawcounts.2$Primer_name, Raw_counts_Test = rawcounts.2$Raw_counts_Test) ###change the order of the columns
rawcounts.2$Primer_name <- gsub(pattern = " ", replacement = "", x = rawcounts.2$Primer_name)
rawcounts.2$Primer_name <- gsub(pattern = "-", replacement = "_", x = rawcounts.2$Primer_name)

plyr::join(rawcounts.1, rawcounts.2, by= "Primer_name")-> rawcounts

rawcounts%>%
  rowwise()%>%
  mutate(Total_reads = sum(Raw_counts_Main, Raw_counts_Test))-> rawcounts

rm(PS.1, PS.1.l, PS.2, rawcounts.1, rawcounts.2, MA.1, MA.2, propMerged.1, propMerged.2, 
   fastqList, filter, filter.track, PipSum, primer, ptable, samplesList, files, filt_path, filtFs,
   filtRs, fullpath, M1, MA, MAList,path, primerF, primerR, samplesAll, fastqFall, fastqRall, trackingF)
