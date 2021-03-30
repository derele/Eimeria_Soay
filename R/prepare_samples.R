ParaData <- read.csv("Data/Soay Eimeria spp.csv", comment.char="%")
DNAData <-  read.csv("Data/Eimeria_DNA_extractions_pilot_Nanodrop_EH.csv")

## some empty cells should be ZERO?!
ParaData[is.na(ParaData)] <- 0

table(ParaData$SAMPLE.REF.%in%DNAData$Tag)

ParaData$SAMPLE.REF[!ParaData$SAMPLE.REF%in%DNAData$Tag]

## One row in the parasite data is duplicated
dupSamp <- ParaData$SAMPLE.REF.[duplicated(ParaData$SAMPLE.REF.)]
ParaData[ParaData$SAMPLE.REF.%in%dupSamp, ]

## and they don't agree :-( (wrote and email to Antonio on 30/04/2021)

## for now
ParaData <- ParaData[!ParaData$SAMPLE.REF.%in%dupSamp, ]

rownames(ParaData) <- ParaData$SAMPLE.REF

pdf("Figures/ParaDataHeat.pdf")
pheatmap(t(ParaData[, !colnames(ParaData)%in%c("TOTAL", "SAMPLE.REF.")]))
dev.off()

### For some samples with DNA data we don't have Eimeria counting
DNAData$Tag[!DNAData$Tag%in%ParaData]

SeqData <-  read.csv("Data/Soay_Sheep_Sample_Indexing.csv")

table(SeqData$DNA_Origin, SeqData$Sample_Origin)

## Hmmm 64 Soay samples, 60 from "oocysts", 4 from Sporocysts
table(SeqData[SeqData$Sample_Origin%in%"Soay_Sheep",  "DNA_Origin"])

## The mouse, Eimeria a positive samples could be interesting later,
SeqData[SeqData$Sample_Origin%in%"Mouse",  ]
SeqData[SeqData$Sample_Origin%in%"Eimeria",  ]

## But for now we focus on the Soay!
SoSeqD <- SeqData[SeqData$Sample_Origin%in%"Soay_Sheep", ]

### 32 samples
length(unique(SoSeqD$Tag))
## for all of which we have PCR duplicates!
table(SoSeqD$Tag)
## amplified on different chips
table(SoSeqD$Tag, SoSeqD$Chip_number)

## for 28 of the sequenced samples we have parasite count data (for 4 not)
table(unique(SoSeqD$Tag)%in%ParaData$SAMPLE.REF.)

## this means all 28 counted samples have been PCRed for sequencing 
table(ParaData$SAMPLE.REF.%in%unique(SoSeqD$Tag))



