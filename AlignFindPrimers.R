library(DECIPHER)
library(ggtree)
library(phangorn)
library(pheatmap)
library(parallel)

## Soay Eimeria species from CRAIG et al. doi:10.1017/S0031182006001144
soay <- c("E. ahsata", "E. bakuensis", "E. crandallis", "E. faurei",
          "E. granulosa", "E. intricata", "E. marsica",
          "E. ovinoidalis", "E. pallida", "E. parva",
          "E. weybridgensis")

### importing data from NCBI: first from a blast of E_ovinoidalis full
### length 18S
BL18S <- readDNAStringSet("./NCBI_data/E_ovinoidalis_18S_BlastRes.fasta")

BL18S_Aln <- AlignSeqs(RNAStringSet(BL18S))
orient <- DistanceMatrix(BL18S_Aln)[1, ] > 0.2
BL18S_Aln[orient] <- reverseComplement(BL18S_Aln[orient])
BL18S_Aln <- AlignSeqs(RemoveGaps(BL18S_Aln))

max(DistanceMatrix(BL18S_Aln))


## then when looking for all the species in Taxonomy and downloading
## all 18S seq
TA18S <- readDNAStringSet("./NCBI_data/E_soay_18S_all.fasta")
## only the loger ones are usable... the shorter are too fragmetned
TA18S_Aln <- AlignSeqs(RNAStringSet(TA18S[width(TA18S)>1000]))

max(DistanceMatrix(TA18S_Aln))

ALN <- AlignProfiles(BL18S_Aln, TA18S_Aln)

max(DistanceMatrix(ALN))

writeXStringSet(ALN, "18S_aln.fasta")

masked_ALN <- MaskAlignment(ALN, includeTerminalGaps = TRUE, maxFractionGaps = 0.4)

ALN_clean <- as(masked_ALN, "RNAStringSet")

ALN_clean <- AdjustAlignment(ALN_clean)

writeXStringSet(ALN_clean, "18S_aln_clean.fasta")

names(ALN_clean) <- gsub(".*(Eimeria )(\\w*).*", "E. \\2", names(ALN_clean))

NJtree <- NJ(DistanceMatrix(ALN_clean))

### root(NJtree, )

is.seq <- sapply(soay, function (x) any(grepl(x, names(ALN_clean))))
table(is.seq)

### Plot a tree
pdf("Figures/18Stree.pdf")
ggtree(NJtree, root.position=25) +
    geom_tiplab(aes(color = label %in% soay) )

dev.off()

### The goal here is to test the existing primers. Get them from the
### latest Hyena project
ptable <- read.csv(file = "/SAN/Victors_playground/Metabarcoding/AA_Hyena/primer_list.csv",
                   sep=",", header=TRUE, stringsAsFactors=FALSE)


Amplicons <- lapply(1:nrow(ptable), function (i) {
    AmplifyDNA(DNAStringSet(c(ptable[i, "Seq_F"], ptable[i, "Seq_R"])),
           RemoveGaps(DNAStringSet(ALN_clean)),
           P=0.001, annealingTemp=60, maxProductSize=700,
           minEfficiency = 0.999)
})


primerMatches <- lapply(Amplicons, function (amp) {
    primerMatches <- lapply(as.character(amp), function (x){
        matched <- grepl(x, RemoveGaps(DNAStringSet(ALN_clean)))
        unique(names(ALN_clean)[matched])
    })
})

Diff_seq_primer <- lapply(primerMatches, function (x) {
    unlist(x[unlist(lapply(x, length))==1])
})

names(Diff_seq_primer) <- apply(ptable[, c("Name_F", "Name_R")], 1, paste, collapse="_")

Diff_seq_primer <- Diff_seq_primer[!unlist(lapply(Diff_seq_primer, is.null))]

X18S_idAble <- unique(unlist(Diff_seq_primer))

X18S_idAble[X18S_idAble %in% soay]

is.seq[is.seq]

## we'd confuse with something in EVERY primer pair: 
### only E. ahsata, but with different Soay Eimeria for other
### primers. Therefore that's the only one.

## so the usable 18S primers are: 
names(Diff_seq_primer)

## and here is how good they are, how many species they distinguish:
lapply(Diff_seq_primer, length)


### PRIMERS TARGETING COI #############################################

COI <- readDNAStringSet("./NCBI_data/Eimeria_Soay_COI.fasta")

COI_Aln <- AlignSeqs(COI)

names(COI_Aln) <- gsub(".*(Eimeria )(\\w*).*", "E\\2", names(COI_Aln))

max(DistanceMatrix(COI_Aln))

writeXStringSet(COI_Aln, "COI_aln.fasta")

writeXStringSet(RemoveGaps(COI_Aln), "COI_clean.fasta")

### read the primers in a format ment for
### https://www.bioinformatics.org/sms2/primer_map.html
### see figure "Figures/Primer_Map_COI.png"

CPr <- readLines("COI_Soay_primerMan.txt")
CPr <- lapply(CPr, strsplit, " ")

PrimerCOIvec <- unlist(lapply(lapply(CPr, "[[", 1), "[", 2))

PrimerCOIvec <- gsub(",", "", PrimerCOIvec)

PrimersCOI <- DNAStringSet(PrimerCOIvec)

names(PrimersCOI) <- unlist(lapply(lapply(CPr, "[[", 1), "[", 1))


library(TmCalculator)

sapply(PrimerCOIvec, Tm_Wallace)
sapply(PrimerCOIvec, Tm_NN)
sapply(PrimerCOIvec, Tm_GC)

writeXStringSet(PrimersCOI, "COI_Soay_Primers.fasta")


### DESIGN PRIMERS TARGETING 5.8S ITS2 and (beginning of) 28S
### #############################


BL28S <- readDNAStringSet("./NCBI_data/CowChicken_Eimeria28S.fasta")

names(BL28S) <- gsub(".*(Eimeria )(\\w*).*", "E_\\2", names(BL28S))

X28S_Aln <- AlignSeqs(RNAStringSet(BL28S))

writeXStringSet(X28S_Aln, "28S_aln.fasta")

masked_28SALN <- MaskAlignment(X28S_Aln, includeTerminalGaps = TRUE,
                               maxFractionGaps = 0.1)

X28S_Aln_clean <- as(masked_28SALN, "RNAStringSet")

X28S_Aln_clean <- AdjustAlignment(X28S_Aln_clean)

writeXStringSet(DNAStringSet(X28S_Aln_clean), "28S_aln_clean.fasta")


X28Pr <- readLines("28S_Soay_primer_map.txt")
X28Pr <- lapply(X28Pr, strsplit, " ")

Primer28vec <- unlist(lapply(lapply(X28Pr, "[[", 1), "[", 2))

Primer28vec <- gsub(",", "", Primer28vec)

Primers28 <- DNAStringSet(Primer28vec)

names(Primers28) <- unlist(lapply(lapply(X28Pr, "[[", 1), "[", 1))

writeXStringSet(Primers28, "28S_Soay_Primers.fasta")


#### DECIPHER PRIMER DESIGN FAILS for COI
## Seqs2DB(COI_Aln, "XStringSet", "/SAN/db/Eimeria_soay_COI_aln.sql",
##         identifier="allID") 

## tiles <- TileSeqs("/SAN/db/Eimeria_soay_COI_aln.sql")


## COIMetaPrimers1 <- DesignPrimers(tiles, annealingTemp=60, minProductSize=300,
##                                  maxProductSize=400, numPrimerSets=3,
##                                  minGroupCoverage = 0.9999, maxDistance=0,
##                                  processors=20)

## COIMetaPrimers2 <- DesignPrimers(tiles, annealingTemp=60, minProductSize=300,
##                                  maxProductSize=400, numPrimerSets=3,
##                                  start=140)

## COIMetaPrimers3 <- DesignPrimers(tiles, annealingTemp=60, minProductSize=300,
##                                  maxProductSize=400, numPrimerSets=3,
##                                  start=290)



### DESIGN PRIMERS TARGETING Apicoplast genomes
