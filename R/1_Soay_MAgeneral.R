##Project: Access Array Soay Sheep
##Aim: Sequencing pre-processing and denoising, merge by primer combination.
##Author: Emanuel Heitlinger mod. by Víctor Hugo Jarquín-Díaz

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

## AData from prepare_samples.R
source("R/prepare_samples.R")

## let's keep track of the names of Eimeria Soay species
soay <- c("E. ahsata", "E. bakuensis", "E. crandallis", "E. faurei",
          "E. granulosa", "E. intricata", "E. marsica",
          "E. ovinoidalis", "E. pallida", "E. parva",
          "E. weybridgensis")

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

MS <- addSampleData(MS, AData)

## This should also be automatic!!!
MS <- MS[, colnames(MS)%in%rownames(getSampleData(MS))]

pdf("Figures/pipeline_heatmapInput.pdf", width=16, height=16)
pheatmap(log10(getCounts(getDerepF(M, dropEmpty=FALSE), what="input")+1))
dev.off()

pheatmap(getCounts(getDerepF(M, dropEmpty=FALSE), what="uniques"))

pheatmap(getCounts(getMergers(M, dropEmpty=FALSE), what="uniques"))

## With the renamed MS this doesn't work! Here for MA
pheatmap(log10(getCounts(getSequenceTable(M, dropEmpty=FALSE),
                         what="uniques")+1))

pdf("Figures/pipeline_heatmapOutput.pdf", width=16, height=16)
pheatmap(log10(getCounts(getSequenceTableNoChime(M, dropEmpty=FALSE),
                         what="uniques")+1))
dev.off()

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
              mean_length = mean(nchar(OTU)),
              n_Samples = n_distinct(Sample[Abundance>0]),
              n_Eim_Samples = n_distinct(Sample[Abundance>0 &
                                                genus %in% "Eimeria"]),
              n_Eim_Soay_Samples = n_distinct(Sample[Abundance>0 &
                                                     genus %in% "Eimeria" &
                                                     Sample_Origin %in% "Soay_Sheep"
                                                     ]),
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
              .) %>%
    arrange(desc(n_Api_ASVs)) ->
    primer.table


ptable$onename <- paste(ptable$Name_F, ptable$Name_R, sep=".")

primer.table <-inner_join(ptable, primer.table, by=c("onename" ="primer"))

colnames(primer.table)


theme_set(theme_bw())

pdf("Figures/sizes.pdf")
ggplot(primer.table, aes(mean_length, tot_counts, label=Gen)) +
    geom_point(aes(color=Gen, size=n_ASVs)) +
    scale_y_log10() + 
    stat_smooth(method="lm", se=FALSE) +
    geom_text(check_overlap = FALSE, nudge_y=0.06) 
dev.off()


pdf("Figures/sizesEASVs.pdf")
ggplot(primer.table, aes(mean_length, tot_counts,
                         label=n_Eim_ASVs)) +
    geom_point(aes(color=Gen, size=n_Eim_ASVs)) +
    scale_y_log10() + 
    stat_smooth(method="lm", se=FALSE) +
    geom_text(check_overlap = FALSE, nudge_y=0.06) 
dev.off()



pdf("Figures/sizeSamples.pdf")
ggplot(primer.table, aes(tot_counts, n_Samples)) +
        geom_point(aes(color=Gen, size=n_Eim_ASVs)) 
dev.off()

pdf("Figures/sizeEimSamples.pdf")
ggplot(primer.table, aes(tot_counts, n_Eim_Samples)) +
        geom_point(aes(color=Gen, size=n_Eim_ASVs)) 
dev.off()

pdf("Figures/sizeEimSoaySamples.pdf")
ggplot(primer.table, aes(tot_counts, n_Eim_Soay_Samples)) +
        geom_point(aes(color=Gen, size=n_Eim_ASVs)) 
dev.off()

pdf("Figures/EASVsEspec.pdf")
ggplot(primer.table, aes(n_Eim_ASVs,
                         Eim_species)) +
    geom_point(aes(color=Gen, size=tot_counts)) +
    stat_smooth(method="lm", se=FALSE) 
dev.off()


### So now let's get all the OTUs for the high diversity Eimeria
### amplicon

primer.table %>% filter(n_Eim_ASVs > 55 &
                        Eim_species > 5) 

## let's decide for the primer above...

filter(Ptab, primer%in%"EimSS_131F_157_F.EimSS_473R_160_R") %>%
    filter(genus%in%"Eimeria" &
           Sample_Origin %in% "Soay_Sheep") %>%
    group_by(OTU) %>% 
    summarize(abundance = sum(Abundance),
              nSamples =  n_distinct(Sample[Abundance>0]),
              species = unique(species)) %>%
    filter(abundance>5 & nSamples > 1) %>% 
    transform(species = gsub("imeria ", ". ", species)) %>%
    transform(origin = "SEQ") %>%
    transform(seqName = 1:nrow(.)) ->
    Eim1stFrame

## now together with the 18S aligment from R/AlignFindPrimers.R
db.seq <- as_tibble(
    cbind(OTU=as.character(RemoveGaps(ALN_clean)),
          species=names(ALN_clean)))

db.seq %>%
    transform(abundance = 1) %>%
    transform(nSamples = 1) %>%
    transform(origin = "NCBInt") %>%
    rbind(Eim1stFrame[,colnames(.)]) %>%
    transform(soay = ifelse(species%in%soay, "Soay", "Other")) %>%
    transform(seqName = make.unique(paste0(species, "_", soay))) ->
    Eim1stFrame


Eim1stFrame %>%
    pull(OTU, seqName)   -> Eim18S1stSeq

Eim18S1stSeq <- RNAStringSet(gsub("T", "U", Eim18S1stSeq))

Eim1stAln <- AlignSeqs(Eim18S1stSeq)

## have a look
DistanceMatrix(Eim1stAln)[1, ]

badSeq <- DistanceMatrix(Eim1stAln)[1, ] > 0.2

Eim1stAln <- Eim1stAln[!badSeq]

NJtree <- NJ(DistanceMatrix(Eim1stAln))

Eim1stFrame %>%
    filter(seqName %in% NJtree$tip.label) %>%
    relocate(seqName) ->
    Eim1stFrame


## tip.groups <- NJtree$tip.label
## names(tip.groups) <- Eim1stFrame$soay
## NJtree <- groupOTU(NJtree, tip.groups)

t <- ggtree(NJtree)

p <- facet_plot(t+xlim_tree(0.2), panel='log10(abundance)',
                data=Eim1stFrame, geom=geom_segment,
                aes(x=0, xend=log10(abundance), y=y, yend=y, color=nSamples), 
                size=3) 


### Plot a tree
pdf("Figures/18SEimSoaytree.pdf", height=10, width=16)
p + geom_tiplab(color=ifelse(Eim1stFrame$origin%in%"SEQ", "salmon", "steelblue")) +
    theme_tree2()
dev.off()


### So now let's get all the OTUs for the deeply sequenced average diversity Eimeria
### amplicon

primer.table %>% filter(Gen %in% "18S" &
                        Eim_species == 3) 

## let's decide for the primer above...

filter(Ptab, primer%in%"wang1141_13_F.Nem_0425_6_3_R") %>%
    filter(genus%in%"Eimeria" &
           Sample_Origin %in% "Soay_Sheep") %>%
    group_by(OTU) %>% 
    summarize(abundance = sum(Abundance),
              nSamples =  n_distinct(Sample[Abundance>0]),
              species = unique(species)) %>%
    filter(abundance>5 & nSamples > 1) %>% 
    transform(species = gsub("imeria ", ". ", species)) %>%
    transform(origin = "SEQ") %>%
    transform(seqName = 1:nrow(.)) ->
    Eim2ndFrame

## now together with the 18S aligment from R/AlignFindPrimers.R
db.seq <- as_tibble(
    cbind(OTU=as.character(RemoveGaps(ALN_clean)),
          species=names(ALN_clean)))

db.seq %>%
    transform(abundance = 1) %>%
    transform(nSamples = 1) %>%
    transform(origin = "NCBInt") %>%
    rbind(Eim2ndFrame[,colnames(.)]) %>%
    transform(soay = ifelse(species%in%soay, "Soay", "Other")) %>%
    transform(seqName = make.unique(paste0(species, "_", soay))) ->
    Eim2ndFrame


Eim2ndFrame %>%
    pull(OTU, seqName)   -> Eim18S2ndSeq

Eim18S2ndSeq <- RNAStringSet(gsub("T", "U", Eim18S2ndSeq))

Eim2ndAln <- AlignSeqs(Eim18S2ndSeq)

## have a look
DistanceMatrix(Eim2ndAln)[1, ]

badSeq <- DistanceMatrix(Eim2ndAln)[1, ] > 0.2

Eim2ndAln[badSeq] <- reverseComplement(Eim2ndAln[badSeq])

Eim2ndAln <- AlignSeqs(RemoveGaps(Eim2ndAln))

NJtree <- NJ(DistanceMatrix(Eim2ndAln))

Eim2ndFrame %>%
    filter(seqName %in% NJtree$tip.label) %>%
    relocate(seqName) ->
    Eim2ndFrame


## tip.groups <- NJtree$tip.label
## names(tip.groups) <- Eim2ndFrame$soay
## NJtree <- groupOTU(NJtree, tip.groups)

t <- ggtree(NJtree)

p <- facet_plot(t+xlim_tree(0.2), panel='log10(abundance)',
                data=Eim2ndFrame, geom=geom_segment,
                aes(x=0, xend=log10(abundance), y=y, yend=y, color=nSamples), 
                size=3) 


### Plot a tree
pdf("Figures/18SEimSoaytree2nd.pdf", height=10, width=16)
p + geom_tiplab(color=ifelse(Eim2ndFrame$origin%in%"SEQ", "salmon", "steelblue")) +
    theme_tree2()
dev.off()


### So now let's get all the OTUs for the medium deeply sequenced average diversity Eimeria
### amplicon

primer.table %>% filter(Gen %in% "18S" &
                        Eim_species == 4) 

## let's decide for the primer above...

filter(Ptab, primer%in%"18S_0067a_deg_5Mod_52_F.NSR399_5Mod_52_R") %>%
    filter(genus%in%"Eimeria" &
           Sample_Origin %in% "Soay_Sheep") %>%
    group_by(OTU) %>% 
    summarize(abundance = sum(Abundance),
              nSamples =  n_distinct(Sample[Abundance>0]),
              species = unique(species)) %>%
    filter(abundance>5 & nSamples > 1) %>% 
    transform(species = gsub("imeria ", ". ", species)) %>%
    transform(origin = "SEQ") %>%
    transform(seqName = 1:nrow(.)) ->
    Eim3rdFrame


### as this has shown itself to be the best primer
filter(Ptab, primer%in%"18S_0067a_deg_5Mod_52_F.NSR399_5Mod_52_R") %>%
    group_by(phylum) %>%
    summarize(abundance = sum(Abundance),
              ASVs = n_distinct(species)) %>%
    write.csv(file="best_18S_phyla.csv")


## now together with the 18S aligment from R/AlignFindPrimers.R
db.seq <- as_tibble(
    cbind(OTU=as.character(RemoveGaps(ALN_clean)),
          species=names(ALN_clean)))

db.seq %>%
    transform(abundance = 1) %>%
    transform(nSamples = 1) %>%
    transform(origin = "NCBInt") %>%
    rbind(Eim3rdFrame[,colnames(.)]) %>%
    transform(soay = ifelse(species%in%soay, "Soay", "Other")) %>%
    transform(seqName = make.unique(paste0(species, "_", soay))) ->
    Eim3rdFrame


Eim3rdFrame %>%
    pull(OTU, seqName)   -> Eim18S3rdSeq

Eim18S3rdSeq <- RNAStringSet(gsub("T", "U", Eim18S3rdSeq))

Eim3rdAln <- AlignSeqs(Eim18S3rdSeq)

## have a look
DistanceMatrix(Eim3rdAln)[1, ]

badSeq <- DistanceMatrix(Eim3rdAln)[1, ] > 0.2

Eim3rdAln[badSeq] <- reverseComplement(Eim3rdAln[badSeq])

Eim3rdAln <- AlignSeqs(RemoveGaps(Eim3rdAln))

NJtree <- NJ(DistanceMatrix(Eim3rdAln))

Eim3rdFrame %>%
    filter(seqName %in% NJtree$tip.label) %>%
    relocate(seqName) ->
    Eim3rdFrame


## tip.groups <- NJtree$tip.label
## names(tip.groups) <- Eim3rdFrame$soay
## NJtree <- groupOTU(NJtree, tip.groups)

t <- ggtree(NJtree)

p <- facet_plot(t+xlim_tree(0.2), panel='log10(abundance)',
                data=Eim3rdFrame, geom=geom_segment,
                aes(x=0, xend=log10(abundance), y=y, yend=y, color=nSamples), 
                size=3) 


### Plot a tree
pdf("Figures/18SEimSoaytree3rd.pdf", height=10, width=16)
p + geom_tiplab(color=ifelse(Eim3rdFrame$origin%in%"SEQ", "salmon", "steelblue")) +
    theme_tree2()
dev.off()


primer.table %>% filter(Gen %in% "18S" &
                        Region %in% "V1-V2") 


## let's decide for the primer above...

filter(Ptab, primer%in%"MarkN_10_F.Euk360_CR_21_R") %>%
    filter(genus%in%"Eimeria" &
           Sample_Origin %in% "Soay_Sheep") %>%
    group_by(OTU) %>% 
    summarize(abundance = sum(Abundance),
              nSamples =  n_distinct(Sample[Abundance>0]),
              species = unique(species)) %>%
    filter(abundance>5 & nSamples > 1) %>% 
    transform(species = gsub("imeria ", ". ", species)) %>%
    transform(origin = "SEQ") %>%
    transform(seqName = 1:nrow(.)) ->
    Eim4thFrame

## now together with the 18S aligment from R/AlignFindPrimers.R
db.seq <- as_tibble(
    cbind(OTU=as.character(RemoveGaps(ALN_clean)),
          species=names(ALN_clean)))

db.seq %>%
    transform(abundance = 1) %>%
    transform(nSamples = 1) %>%
    transform(origin = "NCBInt") %>%
    rbind(Eim4thFrame[,colnames(.)]) %>%
    transform(soay = ifelse(species%in%soay, "Soay", "Other")) %>%
    transform(seqName = make.unique(paste0(species, "_", soay))) ->
    Eim4thFrame


Eim4thFrame %>%
    pull(OTU, seqName)   -> Eim18S4thSeq

Eim18S4thSeq <- RNAStringSet(gsub("T", "U", Eim18S4thSeq))

Eim4thAln <- AlignSeqs(Eim18S4thSeq)

## have a look
DistanceMatrix(Eim4thAln)[1, ]

badSeq <- DistanceMatrix(Eim4thAln)[1, ] > 0.2

Eim4thAln[badSeq] <- reverseComplement(Eim4thAln[badSeq])

Eim4thAln <- AlignSeqs(RemoveGaps(Eim4thAln))

NJtree <- NJ(DistanceMatrix(Eim4thAln))

Eim4thFrame %>%
    filter(seqName %in% NJtree$tip.label) %>%
    relocate(seqName) ->
    Eim4thFrame


## tip.groups <- NJtree$tip.label
## names(tip.groups) <- Eim4thFrame$soay
## NJtree <- groupOTU(NJtree, tip.groups)

t <- ggtree(NJtree)

p <- facet_plot(t+xlim_tree(0.2), panel='log10(abundance)',
                data=Eim4thFrame, geom=geom_segment,
                aes(x=0, xend=log10(abundance), y=y, yend=y, color=nSamples), 
                size=3) 


### Plot a tree
pdf("Figures/18SEimSoaytree4th.pdf", height=10, width=16)
p + geom_tiplab(color=ifelse(Eim4thFrame$origin%in%"SEQ", "salmon", "steelblue")) +
    theme_tree2()
dev.off()



primer.table %>% filter(Gen %in% "COI")


## let's decide for the primer above...

filter(Ptab, primer%in%"modFW11_147_F.modRev10_147_R") %>%
    filter(genus%in%"Eimeria" &
           Sample_Origin %in% "Soay_Sheep") %>%
    group_by(OTU) %>% 
    summarize(abundance = sum(Abundance),
              nSamples =  n_distinct(Sample[Abundance>0]),
              species = unique(species)) %>%
    filter(abundance>5 & nSamples > 1) %>% 
    transform(species = gsub("imeria ", ". ", species)) %>%
    transform(origin = "SEQ") %>%
    transform(seqName = 1:nrow(.)) ->
    Eim5thFrame

## now together with the COI aligment from R/AlignFindPrimers.R
db.seq <- as_tibble(
    cbind(OTU=as.character(RemoveGaps(COI_Aln)),
          species=names(COI_Aln)))

db.seq %>%
    transform(abundance = 1) %>%
    transform(nSamples = 1) %>%
    transform(origin = "NCBInt") %>%
    rbind(Eim5thFrame[,colnames(.)]) %>%
    transform(soay = ifelse(species%in%soay, "Soay", "Other")) %>%
    transform(seqName = make.unique(paste0(species, "_", soay))) ->
    Eim5thFrame


Eim5thFrame %>%
    pull(OTU, seqName)   -> EimCOI5thSeq

EimCOI5thSeq <- DNAStringSet(EimCOI5thSeq)

Eim5thAln <- AlignSeqs(EimCOI5thSeq)

## have a look
DistanceMatrix(Eim5thAln)[1, ]

badSeq <- DistanceMatrix(Eim5thAln)[1, ] > 0.2

NJtree <- NJ(DistanceMatrix(Eim5thAln))

Eim5thFrame %>%
    filter(seqName %in% NJtree$tip.label) %>%
    relocate(seqName) ->
    Eim5thFrame


## tip.groups <- NJtree$tip.label
## names(tip.groups) <- Eim5thFrame$soay
## NJtree <- groupOTU(NJtree, tip.groups)

t <- ggtree(NJtree)

p <- facet_plot(t+xlim_tree(0.2), panel='log10(abundance)',
                data=Eim5thFrame, geom=geom_segment,
                aes(x=0, xend=log10(abundance), y=y, yend=y, color=nSamples), 
                size=3) 


### Plot a tree
pdf("Figures/COIEimSoaytree5th.pdf", height=10, width=16)
p + geom_tiplab(color=ifelse(Eim5thFrame$origin%in%"SEQ", "salmon", "steelblue")) +
    theme_tree2()
dev.off()


filter(Ptab, primer%in%"modFW4_141_F.modRev4_145_R") %>%
    filter(genus%in%"Eimeria" &
           Sample_Origin %in% "Soay_Sheep") %>%
    group_by(OTU) %>% 
    summarize(abundance = sum(Abundance),
              nSamples =  n_distinct(Sample[Abundance>0]),
              species = unique(species)) %>%
    filter(abundance>5 & nSamples > 1) %>% 
    transform(species = gsub("imeria ", ". ", species)) %>%
    transform(origin = "SEQ") %>%
    transform(seqName = 1:nrow(.)) ->
    Eim6thFrame

## now together with the COI aligment from R/AlignFindPrimers.R
db.seq <- as_tibble(
    cbind(OTU=as.character(RemoveGaps(COI_Aln)),
          species=names(COI_Aln)))

db.seq %>%
    transform(abundance = 1) %>%
    transform(nSamples = 1) %>%
    transform(origin = "NCBInt") %>%
    rbind(Eim6thFrame[,colnames(.)]) %>%
    transform(soay = ifelse(species%in%soay, "Soay", "Other")) %>%
    transform(seqName = make.unique(paste0(species, "_", soay))) ->
    Eim6thFrame


Eim6thFrame %>%
    pull(OTU, seqName)   -> EimCOI6thSeq

EimCOI6thSeq <- DNAStringSet(EimCOI6thSeq)

Eim6thAln <- AlignSeqs(EimCOI6thSeq)

## have a look
DistanceMatrix(Eim6thAln)[1, ]

NJtree <- NJ(DistanceMatrix(Eim6thAln))

Eim6thFrame %>%
    filter(seqName %in% NJtree$tip.label) %>%
    relocate(seqName) ->
    Eim6thFrame


## tip.groups <- NJtree$tip.label
## names(tip.groups) <- Eim6thFrame$soay
## NJtree <- groupOTU(NJtree, tip.groups)

t <- ggtree(NJtree)

p <- facet_plot(t+xlim_tree(0.2), panel='log10(abundance)',
                data=Eim6thFrame, geom=geom_segment,
                aes(x=0, xend=log10(abundance), y=y, yend=y, color=nSamples), 
                size=3) 


### Plot a tree
pdf("Figures/COIEimSoaytree6th.pdf", height=10, width=16)
p + geom_tiplab(color=ifelse(Eim6thFrame$origin%in%"SEQ", "salmon", "steelblue")) +
    theme_tree2()
dev.off()

primer.table %>% filter(Gen %in% "23S")


## let's decide for the primer above...

filter(Ptab, primer%in%"Mit23S_731F_165_F.Mit23S_1088R_169_R") %>%
    filter(genus%in%"Eimeria" &
           Sample_Origin %in% "Soay_Sheep") %>%
    group_by(OTU) %>% 
    summarize(abundance = sum(Abundance),
              nSamples =  n_distinct(Sample[Abundance>0]),
              species = unique(species)) %>%
    filter(abundance>5 & nSamples > 1) %>% 
    transform(species = gsub("imeria ", ". ", species)) %>%
    transform(origin = "SEQ") %>%
    transform(seqName = 1:nrow(.)) ->
    Eim6thFrame

## now together with the 23S aligment from R/AlignFindPrimers.R
db.seq <- as_tibble( 
    cbind(OTU=as.character(RemoveGaps(X23S_Aln)),
          species=names(X23S_Aln)))

db.seq %>%
    transform(abundance = 1) %>%
    transform(nSamples = 1) %>%
    transform(origin = "NCBInt") %>%
    rbind(Eim6thFrame[,colnames(.)]) %>%
    transform(soay = ifelse(species%in%soay, "Soay", "Other")) %>%
    transform(seqName = make.unique(paste0(species, "_", soay))) ->
    Eim6thFrame


Eim6thFrame %>%
    pull(OTU, seqName)   -> Eim23S6thSeq

Eim23S6thSeq <- DNAStringSet(Eim23S6thSeq)

Eim6thAln <- AlignSeqs(Eim23S6thSeq)

## have a look
DistanceMatrix(Eim6thAln)[1, ]

badSeq <- DistanceMatrix(Eim6thAln)[1, ] > 0.2

NJtree <- NJ(DistanceMatrix(Eim6thAln))

Eim6thFrame %>%
    filter(seqName %in% NJtree$tip.label) %>%
    relocate(seqName) ->
    Eim6thFrame


## tip.groups <- NJtree$tip.label
## names(tip.groups) <- Eim6thFrame$soay
## NJtree <- groupOTU(NJtree, tip.groups)

t <- ggtree(NJtree)

p <- facet_plot(t+xlim_tree(0.2), panel='log10(abundance)',
                data=Eim6thFrame, geom=geom_segment,
                aes(x=0, xend=log10(abundance), y=y, yend=y, color=nSamples), 
                size=3) 

### Plot a tree
pdf("Figures/23SEimSoaytree6th.pdf", height=10, width=16)
p + geom_tiplab(color=ifelse(Eim6thFrame$origin%in%"SEQ", "salmon", "steelblue")) +
    theme_tree2()
dev.off()


