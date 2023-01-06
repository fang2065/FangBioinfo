setwd("~/Documents/GSDC")

library(oncoPredict)
library(tidyverse)

countdata_predict <- read_tsv("countdata77.txt", comment="#")
countdata_predict <- countdata_predict[!duplicated(countdata_predict$Geneid), ]

testExpr <- countdata_predict %>%
  as.data.frame() %>% 
  column_to_rownames("Geneid") %>% # turn the geneid column into rownames
  as.matrix()

seqdata <- read_tsv("tamoxifen_countdata.txt")
sample_ic50 <- read_tsv("IC50_tamoxifen.txt")

trainingExpr <- seqdata %>%
  as.data.frame() %>% 
  column_to_rownames("GeneID") %>% # turn the geneid column into rownames
  as.matrix()

sample_ic50_training <- sample_ic50 %>%
  as.data.frame() %>% 
  column_to_rownames("CosmicID") %>% # turn the geneid column into rownames
  as.matrix()

calcPhenotype(trainingExprData = trainingExpr,
              trainingPtype = sample_ic50_training,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom='rawData',
              cc=TRUE,
              rsq=TRUE)
