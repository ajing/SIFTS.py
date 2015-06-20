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

combined$

trainIndex <- createResample(combined$FTID, times = 1, list = F)

model_ensemble <- svm(is_disease ~ ., data = subset(combined[trainIndex, ], select = -FTID), type = "C-classification", probability = T)

pre_result <- predict(model_ensemble, combined[-trainIndex, ], probability = T)
pre_result_ensemble <- attr(pre_result, "probabilities")[,2]

test_result <- merge(scores, data.frame(FTID = combined[-trainIndex, ]$FTID, Ensemble = pre_result_ensemble), by = "FTID")

scores_melt <- melt(test_result, id=c("FTID", "VarType"))
scores_melt <- subset(scores_melt, VarType %in% c("Polymorphism", "Disease") & !is.na(value))
scores_melt$VarType <- factor(scores_melt$VarType, levels = c("Polymorphism", "Disease"))
scores_melt <- split(scores_melt, scores_melt$variable)


### Plot with test
p <- rocplot.multiple(scores_melt, groupName = "VarType", predName = "value", title = "ROC Plot", p.value = F)
ggsave(filename = "tmp.pdf")
