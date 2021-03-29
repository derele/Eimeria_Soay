library(ggplot2)
library(tidyverse)

EW <- read.csv("./Eimeria_DNA_extractions_pilot_Nanodrop_EH.csv")

## flot.cols <- which(grepl("\\.flotation", colnames(EW)))

## faecs.cols <- which(grepl("\\.faeces", colnames(EW)))

gather(EW, key="matrix", value="A260.A280",
       c("A260.A280.flotation", "A260.A280.faeces")) %>%
    mutate(matrix=gsub("A260\\.A280\\.", "", matrix)) %>%
    mutate(A260 = ifelse(matrix%in%"flotation", A260.flotation,
                         A260.faeces)) %>%
    mutate(A280 = ifelse(matrix%in%"flotation", A280.flotation,
                         A280.faeces)) %>%
    mutate(A260.A230 = ifelse(matrix%in%"flotation", A260.A230.flotation,
                              A260.A230.faeces)) %>%
    mutate(Nanodrop.ng.ul =  ifelse(matrix%in%"flotation", Nanodrop..ng.ul..flotation,
                                    Nanodrop..ng.ul..faeces)) %>%
    select(-matches( "\\.flotation|\\.faeces")) -> EWL


ggplot(EWL, aes(Nanodrop.ng.ul, A260.A230, color=matrix)) +
    geom_point()



ggplot(EWL, aes(Nanodrop.ng.ul, A260.A280, color=matrix)) +
    geom_point() +
    scale_y_continuous(limits=c(0, 5))


ggplot(subset(EWL, matrix%in%"flotation"), aes(Nanodrop.ng.ul,  Oocysts.counted, color=A260.A280)) +
  geom_point() +
  stat_smooth(method="lm", se=FALSE)


