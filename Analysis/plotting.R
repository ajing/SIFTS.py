#
#
#

polybind <- read.table("./Data/polybind.txt", sep = "\t", header = F, quote = "", na.string = "\\N")
colnames(polybind) = c("PDBID", "chainID", "resnum", "UniProtID", "uniprot_resnum", "VarType", "DiseaseName", "ligName", "distance")

paste(subset(polybind, UniProtID == "P00918" & PDBID == "4MO8" & distance < 4)$resnum, collapse = "+")


length(unique(subset(protein_annotate_withsnp, VarType %in% c("Polymorphism", "Disease"))$UniProtID))
length(unique(subset(protein_annotate_withsnp, VarType %in% c("Polymorphism", "Unclassified"))$UniProtID))
length(unique(subset(protein_annotate_withsnp, VarType %in% c("Unclassified", "Disease"))$UniProtID))

library(dplyr)
# get the proteins with all 3 var type

disease =  245
polymor =  580
unclass =  149
dis_pol =  disease + polymor - length(unique(subset(protein_annotate_withsnp, VarType %in% c("Polymorphism", "Disease"))$UniProtID))
pol_unc =  polymor + unclass - length(unique(subset(protein_annotate_withsnp, VarType %in% c("Polymorphism", "Unclassified"))$UniProtID))
dis_unc =  disease + unclass - length(unique(subset(protein_annotate_withsnp, VarType %in% c("Unclassified", "Disease"))$UniProtID))
dis_unc_poly = 41

install.packages("VennDiagram")
library(VennDiagram)

pdf("vennplot.pdf")
venn.plot <- draw.triple.venn(disease, polymor, unclass, dis_pol, pol_unc, dis_unc, dis_unc_poly, category = c("Disease", "Polymorphism", "Unclassified"), fill = c("red", "blue", "orange"), cex = 2, cat.cex=2, cat.pos = c(-40, 32, 180), cat.col=c("red", "blue", "orange"))
dev.off()


114 / length(unique(subset(protein_annotate_withsnp, VarType %in% c("Disease"))$UniProtID))


interest_preference <- function(p_annotate_bs, interested_par){
  aa_change <- factor(unique(p_annotate_bs[,interested_par]))
  vartype = "Disease"
  vartypes   = c("Disease", "Polymorphism", "Unclassified")
  ratio_amino <- data.frame(merge(aa_change, vartypes, all = TRUE))
  colnames(ratio_amino) <- c(interested_par, "vartypes")
  for (vartype in vartypes) {
    for (each in aa_change){
      each = as.character(each)
      table_res <- table(deal_with_na(p_annotate_bs$VarType == vartype),
                         deal_with_na(as.character(p_annotate_bs[,interested_par]) == each))
      table_res <- apply(table_res, 1:2, as.numeric)
      rownames(table_res) <- c(paste("Not", vartype), vartype)
      colnames(table_res) <- c(paste("Not", each), each)
      
      fish_result = fisher.test(table(data.frame(residueName = (p_annotate_bs[,interested_par] == each), Disease = p_annotate_bs$VarType == vartype)))
      ratio_amino[ratio_amino[,interested_par] == each & ratio_amino$vartypes == vartype,"estimate"] = fish_result$estimate
      ratio_amino[ratio_amino[,interested_par] == each & ratio_amino$vartypes == vartype,"ci_low"] = fish_result$conf.int[1]
      ratio_amino[ratio_amino[,interested_par] == each & ratio_amino$vartypes == vartype,"ci_up"] = fish_result$conf.int[2]
      ratio_amino[ratio_amino[,interested_par] == each & ratio_amino$vartypes == vartype,"pvalue"] = fish_result$p.value
    }
  }
  ratio_amino
}
ratio_amino <- interest_preference(protein_annotate_onlysnp, "uniprot_resnam_3d")

+ theme(axis.text.x = element_text( size=8, angle=30))

all_location <- function(p_annotate_bs, interested_par){
  ratio_amino_all  <- interest_preference(p_annotate_bs, interested_par)
  ratio_amino_core <- interest_preference(subset(p_annotate_bs, location == "Core"), interested_par)
  ratio_amino_surf <- interest_preference(subset(p_annotate_bs, location == "Surface"), interested_par)
  ratio_amino_bs <- interest_preference(subset(p_annotate_bs, location == "Binding Site"), interested_par)
  #ratio_amino <- rbind(cbind(ratio_amino_all, location = "All"),cbind(ratio_amino_core, location = "Core"), cbind(ratio_amino_surf, location = "Surface"), cbind(ratio_amino_bs, location = "Binding Site"))
  
  ratio_amino <- rbind(cbind(ratio_amino_core, location = "Core"), cbind(ratio_amino_surf, location = "Surface"), cbind(ratio_amino_bs, location = "Binding Site"))
  
  ratio_amino_sort <- ratio_amino_all[with(ratio_amino_all, order(-estimate)),]
  ratio_amino[, interested_par] <- factor(ratio_amino[, interested_par], levels = as.character(subset(ratio_amino_sort, vartypes == "Disease")[, interested_par]))
  
  ggplot(subset(ratio_amino, vartypes %in% c("Disease", "Polymorphism")), aes_string(x = interested_par, y = "estimate", ymin = "ci_low", ymax = "ci_up")) + geom_pointrange(aes(col = vartypes), position=position_dodge(width=0.30))  + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("") + scale_y_log10() + facet_grid(location~.)
  ggsave(filename = "tmp.pdf", height=6, width=12)
}

all_location(protein_annotate_onlysnp, "uniprot_resnam_3d")

all_location(protein_annotate_onlysnp, "h_prop_change")
all_location(protein_annotate_onlysnp, "volumn_change")

library(splitstackshape)
concat.split(ratio_amino_core, "h_prop_change", "to")

splitcolumn <- function(ratio_amino, colname){
  Pbefore <- paste(colname, "_1", sep = '')
  Pafter  <- paste(colname, "_2", sep = '')
  ratio_amino[, Pbefore] <- sapply(ratio_amino[, colname], function(X){unlist(strsplit(as.character(X), split = " to ", fixed = T))[1]})
  ratio_amino[, Pafter] <- sapply(ratio_amino[, colname], function(X){unlist(strsplit(as.character(X), split = " to ", fixed = T))[2]})
  ratio_amino
}


ratio_amino_core <- interest_preference(subset(protein_annotate_onlysnp, location == "Core"), "h_prop_change")

all_location_mod <- function(p_annotate_bs, interested_par){
  ratio_amino <- interest_preference(p_annotate_bs, interested_par)
  
  Pbefore <- paste(interested_par, "_1", sep = '')
  Pafter  <- paste(interested_par, "_2", sep = '')

  ratio_amino_sort <- ratio_amino[with(ratio_amino, order(-estimate)),]
  ratio_amino[, interested_par] <- factor(ratio_amino[, interested_par], levels = as.character(subset(ratio_amino_sort, vartypes == "Disease")[, interested_par]))
  
  ratio_amino <- splitcolumn(ratio_amino, interested_par)
  
  ggplot(subset(ratio_amino, vartypes %in% c("Disease", "Polymorphism")), aes_string(x = interested_par, y = "estimate", ymin = "ci_low", ymax = "ci_up")) + geom_pointrange(aes(col = vartypes), position=position_dodge(width=0.50), size = 1)  + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("") + scale_size(range = c(2, 10)) + scale_y_log10() + theme_bw() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + facet_grid(paste(Pafter, "~", Pbefore)) 
  ggsave(filename = "tmp.pdf", height=8, width=8)
}
all_location_mod(subset(protein_annotate_onlysnp, location == "Binding Site"), "h_prop_change")

all_location_mod(subset(protein_annotate_onlysnp, location == "Surface"), "h_prop_change")
all_location_mod(subset(protein_annotate_onlysnp, location == "Core"), "h_prop_change")


all_location_mod(subset(protein_annotate_onlysnp, location == "Binding Site"), "volumn_change")
all_location_mod(subset(protein_annotate_onlysnp, location == "Surface"), "volumn_change")
all_location_mod(subset(protein_annotate_onlysnp, location == "Core"), "volumn_change")


maingraph <- function(p_annotate_bs, interested_par) {
  freq_aa_sort <- freq_aa[with(freq_aa, order(Frequency)),]
  
  ratio_amino_all  <- interest_preference(p_annotate_bs, interested_par)
  ratio_amino_sort <- ratio_amino_all[with(ratio_amino_all, order(-estimate)),]
  print(ratio_amino_sort)
  print(as.character(subset(ratio_amino_sort, vartypes == "Disease")[, "uniprot_resnam_3d"]))
  
  ratio_amino_all[, interested_par] <- factor(ratio_amino_all[, interested_par], levels = as.character(subset(ratio_amino_sort, vartypes == "Disease")[, "uniprot_resnam_3d"]))
  
  freq_loc$AAType <- factor(freq_loc$AAType, levels = as.character(subset(ratio_amino_sort, vartypes == "Disease")[, "uniprot_resnam_3d"]))
  freq_aa$AAType <- factor(freq_aa$AAType, levels = as.character(subset(ratio_amino_sort, vartypes == "Disease")[, "uniprot_resnam_3d"]))
  
  p1 <- ggplot(freq_loc, aes(x = AAType, fill = Location)) + 
    geom_bar(aes(y = Frequency),stat="identity", position="dodge") + 
    scale_y_continuous(labels  = percent) + xlab("") + ylab("Frequency")
  
  p2 <- ggplot(freq_aa, aes(x = AAType, fill = VarType)) + 
    geom_bar(aes(y = Frequency),stat="identity", position="dodge") + 
    scale_y_continuous(labels  = percent) + xlab("") + ylab("Frequency")
  
  p3 <- ggplot(subset(ratio_amino_all, vartypes %in% c("Disease", "Polymorphism")), aes_string(x = interested_par, y = "estimate", ymin = "ci_low", ymax = "ci_up")) + geom_pointrange(aes(col = vartypes), position=position_dodge(width=0.30))  + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("") + scale_y_log10()
  
  
  pdf("tmp.pdf", height=6, width=12)
  
  grid.newpage()
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  
  pushViewport(viewport(layout = grid.layout(3, 1)))
  print(p1, vp = vplayout(1, 1))
  print(p2, vp = vplayout(2, 1))
  print(p3, vp = vplayout(3, 1))
  
  dev.off()
}

maingraph(protein_annotate_onlysnp, "uniprot_resnam_3d")



interest_preference_sig <- function(p_annotate_bs, interested_par){
  aa_change <- factor(unique(p_annotate_bs[,interested_par]))
  vartype = "Disease"
  vartypes   = c("Disease", "Polymorphism", "Unclassified")
  ratio_amino <- data.frame(merge(aa_change, vartypes, all = TRUE))
  colnames(ratio_amino) <- c(interested_par, "vartypes")
  for (vartype in vartypes) {
    for (each in aa_change){
      each = as.character(each)
      table_res <- table(deal_with_na(p_annotate_bs$VarType == vartype),
                         deal_with_na(as.character(p_annotate_bs[,interested_par]) == each))
      table_res <- apply(table_res, 1:2, as.numeric)
      rownames(table_res) <- c(paste("Not", vartype), vartype)
      colnames(table_res) <- c(paste("Not", each), each)
      
      fish_result = fisher.test(table(data.frame(residueName = (p_annotate_bs[,interested_par] == each), Disease = p_annotate_bs$VarType == vartype)))
      ratio_amino[ratio_amino[,interested_par] == each & ratio_amino$vartypes == vartype,"estimate"] = fish_result$estimate
      ratio_amino[ratio_amino[,interested_par] == each & ratio_amino$vartypes == vartype,"ci_low"] = fish_result$conf.int[1]
      ratio_amino[ratio_amino[,interested_par] == each & ratio_amino$vartypes == vartype,"ci_up"] = fish_result$conf.int[2]
      ratio_amino[ratio_amino[,interested_par] == each & ratio_amino$vartypes == vartype,"pvalue"] = fish_result$p.value
    }
  }
  # for aa_change
  ratio_amino <- subset(ratio_amino, pvalue < 0.05 & is.finite(estimate) & vartypes %in% c("Disease", "Polymorphism"))
  #print(ratio_amino)
  #print(as.vector(t(subset(ratio_amino, vartypes == "Disease", select = interested_par))))
  #print(as.vector(t(subset(ratio_amino, vartypes == "Polymorphism", select = interested_par))))
  common_change <- union(as.vector(t(subset(ratio_amino, vartypes == "Disease", select = interested_par))), as.vector(t(subset(ratio_amino, vartypes == "Polymorphism", select = interested_par))))
  print(common_change)
  ratio_amino <- subset(ratio_amino, aa_change %in% common_change)
  ratio_amino
}

all_location_sig <- function(p_annotate_bs, interested_par){
  ratio_amino_core <- interest_preference_sig(subset(p_annotate_bs, location == "Core"), interested_par)
  ratio_amino_surf <- interest_preference_sig(subset(p_annotate_bs, location == "Surface"), interested_par)
  ratio_amino_bs <- interest_preference_sig(subset(p_annotate_bs, location == "Binding Site"), interested_par)
  ratio_amino <- rbind(cbind(ratio_amino_core, location = "Core"), cbind(ratio_amino_surf, location = "Surface"), cbind(ratio_amino_bs, location = "Binding Site"))
  
  ratio_amino_sort <- ratio_amino[with(ratio_amino, order(-estimate)),]
  ratio_amino[, interested_par] <- factor(ratio_amino[, interested_par], levels = as.character(subset(ratio_amino_sort, vartypes == "Disease")[, interested_par]))
  
  print(ratio_amino)
  
  ggplot(subset(ratio_amino, vartypes %in% c("Disease", "Polymorphism")), aes_string(x = interested_par, y = "estimate", ymin = "ci_low", ymax = "ci_up")) + geom_pointrange(aes(col = vartypes), position=position_dodge(width=0.30))  + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("") + scale_y_log10() + facet_grid(location~.)
  ggsave(filename = "tmp.pdf", height=6, width=12)
}

all_location_sig(protein_annotate_onlysnp, "aa_change")


