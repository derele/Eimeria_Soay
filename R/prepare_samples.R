ParaData <- read.csv("Data/Soay Eimeria spp.csv", comment.char="%")

DNAData <-  read.csv("Data/Eimeria_DNA_extractions_pilot_Nanodrop_EH.csv")

table(ParaData$SAMPLE.REF.%in%DNAData$Tag)

ParaData$SAMPLE.REF[!ParaData$SAMPLE.REF%in%DNAData$Tag]

DNAData$Tag[!DNAData$Tag%in%ParaData]

write.csv(sumSample, file="Data/samples.csv")
