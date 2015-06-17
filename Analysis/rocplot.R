rocdata <- function(grp, pred){
  # Produces x and y co-ordinates for ROC curve plot
  # Arguments: grp - labels classifying subject status
  #            pred - values of each observation
  # Output: List with 2 components:
  #         roc = data.frame with x and y co-ordinates of plot
  #         stats = data.frame containing: area under ROC curve, p value, upper and lower 95% confidence interval
  
  grp <- as.factor(grp)
  if (length(pred) != length(grp)) {
    stop("The number of classifiers must match the number of data points")
  } 
  
  if (length(levels(grp)) != 2) {
    stop("There must only be 2 values for the classifier")
  }
  
  cut <- unique(pred)
  tp <- sapply(cut, function(x) length(which(pred > x & grp == levels(grp)[2])))
  fn <- sapply(cut, function(x) length(which(pred < x & grp == levels(grp)[2])))
  fp <- sapply(cut, function(x) length(which(pred > x & grp == levels(grp)[1])))
  tn <- sapply(cut, function(x) length(which(pred < x & grp == levels(grp)[1])))
  tpr <- tp / (tp + fn)
  fpr <- fp / (fp + tn)
  roc = data.frame(x = fpr, y = tpr)
  roc <- roc[order(roc$x, roc$y),]
  
  i <- 2:nrow(roc)
  auc <- (roc$x[i] - roc$x[i - 1]) %*% (roc$y[i] + roc$y[i - 1])/2
  
  pos <- pred[grp == levels(grp)[2]]
  neg <- pred[grp == levels(grp)[1]]
  q1 <- auc/(2-auc)
  q2 <- (2*auc^2)/(1+auc)
  se.auc <- sqrt(((auc * (1 - auc)) + ((length(pos) -1)*(q1 - auc^2)) + ((length(neg) -1)*(q2 - auc^2)))/(length(pos)*length(neg)))
  ci.upper <- auc + (se.auc * 0.96)
  ci.lower <- auc - (se.auc * 0.96)
  
  se.auc.null <- sqrt((1 + length(pos) + length(neg))/(12*length(pos)*length(neg)))
  z <- (auc - 0.5)/se.auc.null
  p <- 2*pnorm(-abs(z))
  
  stats <- data.frame (auc = auc,
                       p.value = p,
                       ci.upper = ci.upper,
                       ci.lower = ci.lower
  )
  
  return (list(roc = roc, stats = stats))
}



rocplot.multiple <- function(test.data.list, groupName = "grp", predName = "res", title = "ROC Plot", p.value = TRUE){
  require(plyr)
  require(ggplot2)
  plotdata <- llply(test.data.list, function(x) with(x, rocdata(grp = eval(parse(text = groupName)), pred = eval(parse(text = predName)))))
  plotdata <- list(roc = ldply(plotdata, function(x) x$roc),
                   stats = ldply(plotdata, function(x) x$stats)
  )
  
  if (p.value == TRUE){
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (P=", signif(p.value, 2), ")", sep=""))
  } else {
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (95%CI ", signif(ci.upper, 2), " - ", signif(ci.lower, 2), ")", sep=""))
  }
  
  mod_level = plotdata$stats$.id[order(-plotdata$stats$auc)]
  plotdata$roc$.id = factor(plotdata$roc$.id, levels = mod_level)
  
  mod_level = paste(names(test.data.list), ": ", annotation, sep = "")[order(-plotdata$stats$auc)]
  label_factor = factor(paste(names(test.data.list), ": ", annotation, sep = ""), levels = mod_level)
  
  p <- ggplot(plotdata$roc, aes(x = x, y = y)) +
    geom_line(aes(colour = .id)) +
    geom_abline (intercept = 0, slope = 1) +
    theme_bw() +
    scale_x_continuous("False Positive Rate (1-Specificity)") +
    scale_y_continuous("True Positive Rate (Sensitivity)") +
    scale_colour_brewer(palette="Set1", labels = mod_level) +
    theme(legend.justification=c(1,0), legend.position=c(1,0), legend.title=element_blank(), legend.key = element_blank())

  return(p)
}

p <- rocplot.multiple(scores_melt, groupName = "VarType", predName = "value", title = "ROC Plot", p.value = F)
ggsave(filename = "tmp.pdf")

p <- rocplot.multiple(list(SVM = data.frame(res = pre_result_svm, grp = protein_annotate_onlysnp$is_disease), Logistic = data.frame(res = pre_result_lm, grp = protein_annotate_onlysnp$is_disease), SIFT = data.frame(res = SIFT_result$SIFT_score, grp = SIFT_result$is_disease)), groupName = "grp", predName = "res", title = "ROC Plot", p.value = F)
ggsave(filename = "tmp.pdf")