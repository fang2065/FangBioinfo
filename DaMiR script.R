library(DaMiRseq)
library(tidyverse)
setwd("~/OneDrive - University of Cambridge/RNA_Seq/BaselineRvsNR")
sampleinfo <- read_tsv("samples727.txt")
seqdata <- read_tsv("counts727.txt", comment="#")
seqdata <- seqdata[!duplicated(seqdata$Geneid), ]
countdata <- seqdata %>%
  as.data.frame() %>% 
  column_to_rownames("Geneid") %>% # turn the geneid column into rownames
  select(sampleinfo$sample) %>% # keep sample columns using sampleinfo$Sample
  as.matrix()
covariatedata <- sampleinfo %>%
  as.data.frame() %>% 
  column_to_rownames("sample") 
covariatedata$class <- as.factor(covariatedata$class)  
covariatedata$batch <- as.factor(covariatedata$batch)# Convert character column to factor
SE<-DaMiR.makeSE(countdata, covariatedata)
data_norm <- DaMiR.normalization(SE, minCounts=10, fSample=0.7,
                                 hyper = "yes", th.cv=3)
data_filt <- DaMiR.sampleFilt(data_norm, th.corr=0.7)
sv <- DaMiR.SV(data_filt)
DaMiR.corrplot(sv, covariatedata, type = "pearson",
               sig.level = 0.01)
data_adjust<-DaMiR.SVadjust(data_filt, sv, n.sv= 4)
DaMiR.Allplot(data_filt, colData(data_filt))
DaMiR.Allplot(data_adjust, colData(data_adjust))
set.seed(12345)
data_clean<-DaMiR.transpose(assay(data_adjust))
df<-colData(data_adjust)
data_reduced <- DaMiR.FSelect(data_clean, covariatedata, th.corr=0.4, type = "pearson")
data_reduced <- DaMiR.FReduct(data_reduced$data)
df.importance <- DaMiR.FSort(data_reduced, covariatedata)
selected_features <- DaMiR.FBest(data_reduced, ranking=df.importance, 
                                 autoselect = "yes")
DaMiR.Clustplot(selected_features$data, df)
Classification_res <- DaMiR.EnsembleLearning(selected_features$data,
                                             classes=covariatedata$class, fSample.tr = 0.7,
                                             fSample.tr.w = 0.7, iter = 100)
#Dataset for prediction 
set.seed(10101)

# Variables initialization
cv_models <- list()
cv_predictors <- list()
res_df <- data.frame(matrix(nrow = nfold, ncol =7))
colnames(res_df) <- c("Accuracy",
                      "N.predictors",
                      "MCC",
                      "sensitivity", 
                      "specificity",
                      "PPV",
                      "NPV")

For i in 
{
nSampl_cl1 <- 4
nSampl_cl2 <- 8

idx_test_cl1 <- sample(1:5, nSampl_cl1)
idx_test_cl2 <- sample(1:10, nSampl_cl2) + 5
idx_test <- c(idx_test_cl1, idx_test_cl2)
Test_set <- data_adjust[, -idx_test, drop=FALSE]
Learning_set <- data_adjust[, idx_test, drop=FALSE]
# Training and Test into a '10fold' Cross Validation
nfold <- 4
cv_sample <- c(rep(seq_len(nfold), each = ncol(Learning_set)/(2*nfold)),
               rep(seq_len(nfold), each = ncol(Learning_set)/(2*nfold)))


  for (cv_fold in seq_len(nfold)) {
    # create Training and validation SET2
    idx_cv <- which(cv_sample != cv_fold)
    TR_set <- Learning_set[,idx_cv, drop = FALSE]
    Val_set <- Learning_set[,-idx_cv, drop = FALSE]
  
    #Feature selection 
    TR_set_clean <- DaMiR.transpose(assay(TR_set))
    df<-colData(TR_set)
    data_reduced <- DaMiR.FSelect(TR_set_clean, as.data.frame(colData(TR_set)), th.corr=0.6)
    data_reduced <- DaMiR.FReduct(data_reduced$data, th.corr = 0.9)
    df.importance <- DaMiR.FSort(data_reduced,
                               as.data.frame(colData(TR_set)))

    selected_features <- DaMiR.FBest(data_reduced,
                                   ranking = df.importance,
                                   n.pred = 5)
  
    #update datasets
    TR_set <- TR_set[selected_features$predictors, drop = FALSE]
    Val_set <- Val_set[selected_features$predictors, drop = FALSE]
  
    #Model Building 
    ensl_model <- DaMiR.EnsL_Train(TR_set)
    cv_models[[cv_fold]] <- ensl_model
  
    #Model testing
    res_Val <- DaMiR.EnsL_Test(Val_set,
                             EnsL_model = ensl_model)
    # Store all ML results
    res_df[cv_fold,1] <- res_Val$accuracy[1]
    res_df[cv_fold,2] <- length(res_Val$predictors)
    res_df[cv_fold,3] <- res_Val$MCC[1]
    res_df[cv_fold,4] <- res_Val$sensitivity[1]
    res_df[cv_fold,5] <- res_Val$Specificty[1]
    res_df[cv_fold,6] <- res_Val$PPV[1]
    res_df[cv_fold,7] <- res_Val$NPV[1]
  
    cv_predictors[[cv_fold]] <- res_Val$predictors
  }
  
  }
