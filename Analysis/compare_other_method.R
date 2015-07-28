#### How is it compared to other methods?

scores <- read.table("/users/ajing/SNPDist/Data/dscore3.1.txt", sep = "\t", header = F, quote = "", na.string = "\\N")

colnames(scores) <- c("FTID","VarType","SIFT","Polyphen_HDIV","Polyphen2","LRT","MutationTaster","MutationAssessor","FATHMM","PROVEAN","MetaSVM")


#astring = ""
#for(i in 1:6){
#  astring = paste(astring, ", 'score", i, "'", ", 'pred", i, "'", sep = "")
#}

library(reshape2)
#scores_melt <- melt(subset(scores[1:10,], select = - grep("pred",colnames(scores))), id=c("FTID", "VarType"))

scores$Polyphen_HDIV <- NULL
scores$MetaSVM <- NULL
scores$MutationTaster <- NULL

scores_melt <- melt(scores[1:100,], id=c("FTID", "VarType"))
scores_melt <- subset(scores_melt, VarType %in% c("Polymorphism", "Disease") & !is.na(value))
scores_melt$VarType <- factor(scores_melt$VarType, levels = c("Polymorphism", "Disease"))
scores_melt <- split(scores_melt, scores_melt$variable)

p <- rocplot.multiple(scores_melt, groupName = "VarType", predName = "value", title = "ROC Plot", p.value = F)
ggsave(filename = "tmp.pdf")


# how does it work if combine with one of other program
library(caret)
trainIndex <- createResample(protein_annotate_onlysnp$VarType, times = 1, list = F)

### Train the model
######## Linear Model

model <- lm(is_disease ~ hydro_change + location + location:size_change, data = protein_annotate_onlysnp[trainIndex, ])

pre_result_lm <- predict(model, protein_annotate_onlysnp[-trainIndex, ])

######## SVM
library('e1071')
model <- svm(is_disease ~ ., data = subset(protein_annotate_onlysnp[trainIndex, ], select = c("is_disease", "hydro_change", "size_change", "location")), type = "C-classification", probability = T)
pre_result <- predict(model, subset(protein_annotate_onlysnp[-trainIndex, ], select = c( "hydro_change", "size_change", "location")), probability = T)
pre_result_svm <- attr(pre_result, "probabilities")[,2]

######### Ensemble
pre_result <- predict(model, subset(protein_annotate_onlysnp[trainIndex, ], select = c( "hydro_change", "size_change", "location")), probability = T)
pre_result <- attr(pre_result, "probabilities")[,2]
combined = merge(subset(scores, select = c("FTID", "FATHMM")), data.frame(FTID = protein_annotate_onlysnp[trainIndex, ]$FTID, is_disease = protein_annotate_onlysnp[trainIndex, ]$is_disease, SVM = pre_result), by = "FTID")
model_ensemble <- svm(is_disease ~ ., data = subset(combined, select = -FTID), type = "C-classification", probability = T)


combined_test = merge(subset(scores, select = c("FTID", "FATHMM")), data.frame(FTID = protein_annotate_onlysnp[-trainIndex, ]$FTID, is_disease = protein_annotate_onlysnp[-trainIndex, ]$is_disease, SVM = pre_result_svm), by = "FTID")
pre_result <- predict(model_ensemble, subset(combined_test, !is.na(FATHMM)), probability = T)
pre_result_ensemble <- attr(pre_result, "probabilities")[,2]

### Test the model
test_result <- scores
#test_result <- merge(test_result, data.frame(FTID = protein_annotate_onlysnp[-trainIndex, ]$FTID, LinearModel = pre_result_lm), by = "FTID")
test_result <- merge(test_result, data.frame(FTID = subset(combined_test, !is.na(FATHMM))$FTID, SVM = subset(combined_test, !is.na(FATHMM))$SVM, Ensemble = 1- pre_result_ensemble), by = "FTID")

scores_melt <- melt(test_result, id=c("FTID", "VarType"))
scores_melt <- subset(scores_melt, VarType %in% c("Polymorphism", "Disease") & !is.na(value))
scores_melt$VarType <- factor(scores_melt$VarType, levels = c("Polymorphism", "Disease"))
scores_melt <- split(scores_melt, scores_melt$variable)

## the result is not so good

############ another way to ensemble
combined = merge(subset(scores, !is.na(FATHMM), select = c("FTID", "FATHMM")), subset(protein_annotate_onlysnp, select = c("FTID", "is_disease", "hydro_change", "size_change", "location")), by = "FTID")

trainIndex <- createResample(combined$FTID, times = 1, list = F)

model_ensemble <- svm(is_disease ~ ., data = subset(combined[trainIndex, ], select = -FTID), type = "C-classification", probability = T)

pre_result <- predict(model_ensemble, combined[-trainIndex, ], probability = T)
pre_result_ensemble <- attr(pre_result, "probabilities")[,2]

test_result <- merge(scores, data.frame(FTID = combined[-trainIndex, ]$FTID, Ensemble = pre_result_ensemble), by = "FTID")

scores_melt <- melt(test_result, id=c("FTID", "VarType"))
scores_melt <- subset(scores_melt, VarType %in% c("Polymorphism", "Disease") & !is.na(value))
scores_melt$VarType <- factor(scores_melt$VarType, levels = c("Polymorphism", "Disease"))
scores_melt <- split(scores_melt, scores_melt$variable)

## still not so good (AUC: 0.85, FATHMM is 0.9)

############# OK, I will use xgboost to kick your ass
# Convert from classes to numbers
require(xgboost)
require(methods)
require(data.table)
require(magrittr)

combined = merge(subset(scores, !is.na(FATHMM), select = c("FTID", "FATHMM")), subset(protein_annotate_onlysnp, select = c("FTID", "VarType", "hydro_change", "size_change", "location")), by = "FTID")
combined = subset(combined, VarType %in% c("Disease", "Polymorphism"))

combined$is_disease[combined$VarType == "Disease"] = 1
combined$is_disease[combined$VarType == "Polymorphism"] = 0
combined$VarType = NULL

combined_xg = subset(combined, select = -FTID)
combined_xg$location[combined_xg$location == "Core"] = 1
combined_xg$location[combined_xg$location == "Surface"] = 2
combined_xg$location[combined_xg$location == "Binding Site"] = 3
combined_xg$location = as.numeric(combined_xg$location)

trainIndex <- createResample(combined_xg$is_disease, times = 1, list = F)

train = data.table(combined_xg[trainIndex, ])
test  = data.table(combined_xg[-trainIndex, ])
y = train$is_disease
train$is_disease = NULL

trainMatrix <- train[,lapply(.SD,as.numeric)] %>% as.matrix
testMatrix <- test[,lapply(.SD,as.numeric)] %>% as.matrix

numberOfClasses <- max(y) + 1

param <- list("objective" = "binary:logistic",
              "eval_metric" = "auc")

cv.nround <- 1000
cv.nfold <- 3

bst.cv = xgb.cv(param=param, data = trainMatrix, label = y, 
                nfold = cv.nfold, nrounds = cv.nround)

# train the real model
nround = which(bst.cv$test.auc.mean == max(bst.cv$test.auc.mean))
# 688
bst = xgboost(param=param, data = trainMatrix, label = y, nrounds=nround)

pre_result_ensemble <- predict(bst, testMatrix)

test_result <- merge(subset(scores, select = -VarType), data.frame(FTID = combined[-trainIndex, ]$FTID, Ensemble = pre_result_ensemble, VarType = combined[-trainIndex, ]$is_disease), by = "FTID")

######## SVM for only three variables
library('e1071')
library(doMC)
registerDoMC(cores = 4)

model <- svm(is_disease ~ hydro_change + location + location:size_change, data = cbind(trainMatrix, is_disease = y), type = "C-classification", probability = T)
pre_result <- predict(model, testMatrix, probability = T)
pre_result_svm <- attr(pre_result, "probabilities")[,2]
# 0.665

fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  classProbs = T,
  repeats = 10)

svmFit <- train(is_disease ~ ., data = data.frame(trainMatrix, is_disease = c("Polymorphism", "Disease")[y + 1]),
                method = "svmLinear",
                trControl = fitControl,
                preProc = c("center", "scale"),
                tuneLength = 8,
                verbose = T,
                metric = "ROC")

pre_result <- predict(model, subset(testMatrix, select = -is_disease), probability = T)
pre_result_svm <- attr(pre_result, "probabilities")[,2]

# adding svm
test_result <- merge(subset(scores, select = -VarType), data.frame(FTID = combined[-trainIndex, ]$FTID, SVMSimple = 1- pre_result_svm, Ensemble = pre_result_ensemble, VarType = combined[-trainIndex, ]$is_disease), by = "FTID")


scores_melt <- melt(test_result, id=c("FTID", "VarType"))
scores_melt <- subset(scores_melt, !is.na(value))
#scores_melt <- subset(scores_melt, VarType %in% c("Polymorphism", "Disease") & !is.na(value))
#scores_melt$VarType <- factor(scores_melt$VarType, levels = c("Polymorphism", "Disease"))
scores_melt <- split(scores_melt, scores_melt$variable)

### Plot with test
p <- rocplot.multiple(scores_melt, groupName = "VarType", predName = "value", title = "ROC Plot", p.value = F)
ggsave(filename = "tmp.pdf")
