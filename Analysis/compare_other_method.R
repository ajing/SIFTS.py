#### How is it compared to other methods?

scores <- read.table("/users/ajing/SNPDist/Data/dscore2.1.txt", sep = "\t", header = F, quote = "", na.string = "\\N")

colnames(scores) <- c("FTID", "VarType", 'score1', 'pred1', 'score2', 'pred2', 'score3', 'pred3', 'score4', 'pred4', 'score5', 'pred5', 'score6', 'pred6')


#astring = ""
#for(i in 1:6){
#  astring = paste(astring, ", 'score", i, "'", ", 'pred", i, "'", sep = "")
#}

library(reshape2)
scores_melt <- melt(subset(scores[1:100,], select = - grep("pred",colnames(scores))), id=c("FTID", "VarType"))
scores_melt <- subset(scores_melt, VarType %in% c("Polymorphism", "Disease") & !is.na(value))
scores_melt$VarType <- factor(scores_melt$VarType, levels = c("Polymorphism", "Disease"))
scores_melt <- split(scores_melt, scores_melt$variable)

p <- rocplot.multiple(scores_melt, groupName = "VarType", predName = "value", title = "ROC Plot", p.value = F)
ggsave(filename = "tmp.pdf")


# how does it work if combine with one of other program
library(caret)
trainIndex <- createResample(scores$VarType,list = F)

### Train the model
######## Linear Model

model <- lm(is_disease ~ hydro_change + location + location:size_change, data = subset(protein_annotate_onlysnp, FTID %in% scores$FTID[trainIndex]))

pre_result_lm <- predict(model, data = subset(protein_annotate_onlysnp, FTID %in% scores$FTID[-trainIndex]))

######## SVM
library('e1071')
model <- svm(is_disease ~ ., data = subset(subset(protein_annotate_onlysnp, FTID %in% scores$FTID[trainIndex]), select = c("is_disease", "hydro_change", "size_change", "location")), type = "C-classification", probability = T)
pre_result <- predict(model, data = subset(subset(protein_annotate_onlysnp, FTID %in% scores$FTID[-trainIndex]), select = c( "hydro_change", "size_change", "location")), probability = T)
pre_result_svm <- attr(pre_result, "probabilities")[,2]

### Test the model
test_result <- subset(scores[-trainIndex,], select = - grep("pred",colnames(scores)))
test_result <- merge(test_result, data.frame(FTID = subset(protein_annotate_onlysnp, FTID %in% scores$FTID[-trainIndex])$FTID, pre_result_lm), by = "FTID", all.x = T)
test_result <- merge(test_result, data.frame(FTID = subset(protein_annotate_onlysnp, FTID %in% scores$FTID[-trainIndex])$FTID, pre_result_svm), by = "FTID", all.x = T)

scores_melt <- melt(test_result, id=c("FTID", "VarType"))


### Plot with test
scores_melt
