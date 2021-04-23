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

doMultiAmpSort <- FALSE

doMultiAmpError <- FALSE

doMultiAmpPipe <- FALSE

doTax <- FALSE

doPS <- FALSE

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


## ## some very interesting patterns in the 
## pheatmap(getCounts(getDerepF(MA, dropEmpty=FALSE), what="input"))

## pheatmap(getCounts(getDerepF(MA, dropEmpty=FALSE), what="uniques"))

## pheatmap(getCounts(getMergers(MA, dropEmpty=FALSE),
## what="uniques"))

## pheatmap(getCounts(getSequenceTableNoChime(MA, dropEmpty=FALSE),
## what="uniques"))

## pheatmap(log10(getCounts(getSequenceTableNoChime(MA,
## dropEmpty=FALSE), what="uniques")+0.1))


##Lets run BLAST
if (doTax) {
  M <- blastTaxAnnot(MA,
                     ### against a very old database
                     db = "/SAN/db/blastdb/nt/nt",
                     negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                     infasta =
                         "/SAN/Victors_playground/AA_Soay/output/Soay_1_in.fasta",
                     outblast =
                         "/SAN/Victors_playground/AA_Soay/output/Soay_1_out.blt",
                     taxonSQL = "/SAN/db/taxonomy/taxonomizr.sql", 
                     num_threads = 20)
  saveRDS(M, file="/SAN/Victors_playground/AA_Soay/tmp/Mtaxed.Rds")
}else{
    M <- readRDS(file="/SAN/Victors_playground/AA_Soay/tmp/Mtaxed.Rds")
}


## the second line should automatically adjust in the whole data
## structure in the package?!!
MS <- M
rownames(MS@sampleData) <- getSampleData(MS)$sampleID
colnames(MS) <- rownames(getSampleData(MS))

## AData from prepare_samples.R <- SOURCE ME!!!
source("R/prepare_samples.R")

MS <- addSampleData(MS, AData)

## This should also be automatic!!!
MS <- MS[, colnames(MS)%in%rownames(getSampleData(MS))]

pheatmap(getCounts(getDerepF(M, dropEmpty=FALSE), what="input"))

pheatmap(getCounts(getDerepF(M, dropEmpty=FALSE), what="uniques"))

pheatmap(getCounts(getMergers(M, dropEmpty=FALSE), what="uniques"))

## With the renamed MS this doesn't work! Here for MA
pheatmap(log10(getCounts(getSequenceTable(M, dropEmpty=FALSE),
                         what="uniques")+1))

pheatmap(log10(getCounts(getSequenceTableNoChime(M, dropEmpty=FALSE),
                         what="uniques")+1))

## ## Working with M for now, which is missing the additional sample Data
##   plotAmpliconNumbers(M)


if (doPS) {
    ##Create individual Phyloseq objects 
    PS <- toPhyloseq(M, samples=colnames(M))
    ##Separated list of primers
    PS.l <- toPhyloseq(M, samples=colnames(M), multi2Single = FALSE) 


    addSD <- function (ps) {
        S <- sample_data(ps)
        attr(S, "class") <- "data.frame"
        SM <- merge(S, AData, by.x="sampleID", by.y="row.names", all.x=TRUE)
        rownames(SM) <- SM$sampleID
        O <- otu_table(ps)
        rownames(O) <- gsub("2021_16_soay_Main_Run_", "", rownames(O))
        SM <- SM[rownames(O),]
        tt <- slot(ps, "tax_table")
        if(!is.null(tt) & all(dim(tt))>0){
            phyloseq(otu_table(O),
                     sample_data(SM),
                     tax_table(tt))
        } else {
            phyloseq(otu_table(O),
                     sample_data(SM))
        }
    }

    PS <- addSD(PS)

    P.l <- lapply(PS.l, function (x) {    
        S <- x
        if(!is.null(S)&length(S)>0){
            addSD(S)
        } else {S}
    })

    P.l <- P.l[!unlist(lapply(P.l, is.null))]

    ## melt the phyloseq object for each amplicon
    P.m <- mclapply(P.l, function (x){
        if(!is.null(x)){
            phyloseq::psmelt(x)
        } else {NA}
    }, mc.cores=10)
    ## saves for the RAM problem on harriet
    saveRDS(PS, file="/home/ele/PS.Rds")
    saveRDS(P.m, file="/home/ele/Pm.Rds")
    saveRDS(M, file="/home/ele/M.Rds")
} else {
    PS <- readRDS(file="/home/ele/PS.Rds")
    P.m <- readRDS(file="/home/ele/Pm.Rds")
    M  <- readRDS(M, file="/home/ele/M.Rds")
}

        
        
### investigate how in all world it can happen that only some
### annotation levels are missing!
annTaxa <- c("superkingdom", "phylum", "order", "class", "family", "genus", "species")

## dumping the stuff without tax annotation for now
has.tax.an <- unlist(lapply(P.m, function (x) all(annTaxa%in%colnames(x))))
names(P.m)[!has.tax.an]
## the ORF470 and one Eimeria 28S primer didn't get any tax annotation
## but also very few reads
lapply(P.m[!has.tax.an], function(x) sum(x$Abundance))
lapply(P.m[!has.tax.an], function(x) x[x$Abundance>0, "Sample"])
##  ORF primers amplified only in samples named PB02, PB01, PB04, PB05
## and the E. falcifomris controls

##  28S primer amplified in a few more of the Soay and the
## E. falcifomris controls: on further inspection 
## mod28F1_149_F.mod28R5_153_R:
table(P.m[["mod28F1_149_F.mod28R5_153_R"]][, "species"])
## has only weird fungus amplified

## Add the primer name to each molten df
Pn <- lapply(seq_along(P.m), function (i){
    cbind(P.m[[i]], primer=names(P.m)[[i]])
})

## and remove the ones without tax annotation (the ORF470 ones)
Pn <- Pn[has.tax.an]

### tidyverse FUN!
library(tidyverse)
## put all in one giant tibble
Ptab <- as_tibble(Reduce(rbind, Pn))

Ptab %>%
    group_by(primer) %>%
    summarise(n_ASVs = n_distinct(OTU),
              n_Samples = n_distinct(Sample[Abundance>0]),
              tot_counts = sum(Abundance),
              tot_species = n_distinct(species),
              tot_genera = n_distinct(genus),
              tot_families = n_distinct(family),
              tot_orders = n_distinct(order),
              tot_class = n_distinct(class),
              tot_phyla = n_distinct(phylum),
              n_Api_ASVs = n_distinct(OTU[phylum%in%"Apicomplexa"]),
              n_Api_Samples = n_distinct(Sample[Abundance>0 &
                                                phylum%in%"Apicomplexa"]),
              Api_counts = sum(Abundance[phylum%in%"Apicomplexa"]),
              Api_species = n_distinct(species[phylum%in%"Apicomplexa"]),
              Api_prop = sum(Abundance[phylum%in%"Apicomplexa"])/sum(Abundance),
              Eim_species = n_distinct(species[genus%in%"Eimeria"]),
              n_Eim_ASVs = n_distinct(OTU[genus%in%"Eimeria"]),
              ) %>%
    arrange(desc(n_Api_ASVs)) ->
    primer.table

colnames(primer.table)

ggplot(primer.table, aes(n_Eim_ASVs, n_Samples,
                         size=log(tot_counts), color=tot_phyla)) +
    geom_point()



countASVsby <- function (x, what){
        group_by(!!what, primer) %>%
        summarise(n_ASVs = n_distinct(OTU),
                  n_Samples = n_distinct(Sample[Abundance>0]),
                  tot_counts = sum(Abundance)) %>%
        arrange(desc(n_ASVs)) 
}


countASVsbyList <- function (L, what){
    LL <- lapply(L, function (x) {
        if(what%in%colnames(x)){
            countASVsby(x, sym(what))
        } else{NA}
    })
    LLg <- LL[unlist(lapply(LL, is.data.frame))]
    ## add the name to the df
    LLg <- lapply(seq_along(LLg), function (i){
        cbind(LLg[[i]], primer=names(LLg)[[i]])
    })
    Reduce(rbind, LLg)
}

    

phylum.tab <- countASVsbyList(P.m, "phylum")
family.tab <- countASVsbyList(P.m, "family")
species.tab <- countASVsbyList(P.m, "species")


api.primer <- subset(phylum.tab, phylum%in%"Apicomplexa")[, "primer"]
nem.primer <- subset(phylum.tab, phylum%in%"Nematoda")[, "primer"]

## ggplot(subset(family.tab, primer%in%api.primer), 

as_tibble(phylum.tab) %>%
    group_by(primer) %>%
    summarise(Api_ASVs = sum(n_ASVs[phylum%in%"Apicomplexa"]),
              n_Samples = n_distinct(Sample[Abundance>0]),
              tot_counts = sum(Abundance)) %>%
       


soay <- c("E. ahsata", "E. bakuensis", "E. crandallis", "E. faurei",
          "E. granulosa", "E. intricata", "E. marsica",
          "E. ovinoidalis", "E. pallida", "E. parva",
          "E. weybridgensis")

soayE <- gsub("E\\.", "Eimeria", soay)

as_tibble(Eimeriadf) %>%
    group_by(species) %>% summarize(Sabu=sum(Abundance)) %>%
    arrange(Sabu) %>% transform(isSheep=species%in%soayE) %>% as.data.frame()
