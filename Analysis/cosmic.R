# how about cancer mutation? does it also likely to be in protein ligand binding sites?

#################################################
#### For cosmic
#################################################
cosmic <- read.table("/users/ajing/SNPDist/Data/CosmicCompleteExport.tsv", sep = "\t", header = T, quote = "", na.string = "")

# Mutation in coding region
library("stringr")
cosmic_missense <- subset(cosmic, Mutation.Description == "Substitution - Missense")
cosmic_missense$Mutation.AA.Before <- apply(cosmic_missense, 1, function(x){str_sub(x[["Mutation.AA"]], 3, 3)}) 
cosmic_missense$Mutation.AA.After <- apply(cosmic_missense, 1, function(x){str_sub(x[["Mutation.AA"]], -1, -1)}) 
cosmic_missense$Mutation.AA.resnum <- apply(cosmic_missense, 1, function(x){as.numeric(str_sub(x[["Mutation.AA"]], 4, -2))}) 

aa_list <- c("R", "Q", "K", "V", "C", "Y", "D", "P", "I", "L", "A", "N", "H", "M", "T", "F", "W", "S", "E", "G", "Z")
cosmic_missense <- subset(cosmic_missense, !is.na(Mutation.AA.resnum) & Mutation.AA.Before %in% aa_list & Mutation.AA.After %in% aa_list)

cosmic_missense$Mutation.AA.Before <- factor(cosmic_missense$Mutation.AA.Before)
cosmic_missense$Mutation.AA.After <- factor(cosmic_missense$Mutation.AA.After)


library("sqldf")
protein_annotate_cancer <- sqldf(c("create index idx on cosmic_missense (`Gene.name`, `Mutation.AA.resnum`)", "select * from protein_annotate left join cosmic_missense on `Gene.name` = protein_annotate.gene_name_acc and `Mutation.AA.resnum` = protein_annotate.resnum"))
protein_annotate_cancer$VarType = NA
protein_annotate_cancer[!is.na(protein_annotate_cancer$Mutation.AA.Before), "VarType"] = "Cancer"

table(protein_annotate_cancer$VarType, protein_annotate_cancer$location)

# more columns for cancer mutations
require(dplyr)
tmp = protein_annotate_withsnp %>% select(uniprot_resnam, uniprot_resnam_3d) %>% distinct  # translate from 1digit amino acide code to 3 digit.

protein_annotate_withcancer = subset(protein_annotate_cancer, UniProtID %in% factor(unique(subset(protein_annotate_cancer, !is.na(VarType))$UniProtID)))

# take unique SNP
library(sqldf)
protein_annotate_withcancer_u = sqldf("select * from protein_annotate_withcancer group by UniProtID, uniprot_resnum, `Mutation.ID`, `Primary.histology`")

protein_annotate_cancer_t = merge(protein_annotate_withcancer_u, tmp, by.x = "Mutation.AA.After", by.y = "uniprot_resnam", suffixes = c(".before",".after"), all.x = T)


# The first table
get_stat_eachtype <- function(p_annotate, snp_type){
  require(dplyr)
  tmp = p_annotate %>% select(UniProtID, uniprot_resnum, location) %>% distinct
  allres_table <- table(tmp$location)
  loc_table = table(subset(p_annotate, VarType == snp_type)$location)
  print("Observed")
  print(loc_table)
  print(sum(loc_table))
  print("Expected:")
  exp_table = allres_table / sum(allres_table) * sum(loc_table)
  print(exp_table)
  print("O/E")
  print(loc_table/exp_table)
}
get_stat_eachtype(protein_annotate_withcancer, "Cancer")

get_stat_eachtype(protein_annotate_withcancer_u, "Cancer")


odds_ratio_stat <- function(p_annotate, vartype){
  require("epitools")
  table_res <- table(deal_with_na(p_annotate$VarType == vartype), deal_with_na(p_annotate$location == "Core"))
  table_res <- apply(table_res, 1:2, as.numeric)
  rownames(table_res) <- c(paste("Not", vartype), vartype)
  colnames(table_res) <- c("Not Core", "Core")
  print(oddsratio.wald(table_res))
  
  table_res <- table(deal_with_na(subset(p_annotate, location %in% c("Surface", "Binding Site"))$VarType == vartype),
                     deal_with_na(subset(p_annotate, location %in% c("Surface", "Binding Site"))$location == "Surface"))
  table_res <- apply(table_res, 1:2, as.numeric)
  rownames(table_res) <- c(paste("Not", vartype), vartype)
  colnames(table_res) <- c("Binding Site", "Surface")
  print(oddsratio.wald(table_res))
  
  table_res <- table(deal_with_na(subset(p_annotate, location %in% c("Core", "Binding Site"))$VarType == vartype),
                     deal_with_na(subset(p_annotate, location %in% c("Core", "Binding Site"))$location == "Core"))
  table_res <- apply(table_res, 1:2, as.numeric)
  rownames(table_res) <- c(paste("Not", vartype), vartype)
  colnames(table_res) <- c("Binding Site", "Core")
  print(oddsratio.wald(table_res))
}

# Cancer
odds_ratio_stat(protein_annotate_withcancer, "Cancer")

odds_ratio_stat(protein_annotate_withcancer_u, "Cancer")


# frequency distribution
freq_plot <- function(protein_annotate_withsnp) {
  require(ggplot2)
  require(reshape2)
  require(scales)
  require(dplyr)
  tmp = protein_annotate_withsnp %>% select(UniProtID, uniprot_resnum, uniprot_resnam_3d.before) %>% distinct
  allres_table <- table(tmp$uniprot_resnam_3d.before)
  
  freq_aa = melt(cbind(prop.table(table(protein_annotate_withsnp$uniprot_resnam_3d.before, protein_annotate_withsnp$VarType), margin = 2),All_Residues = prop.table(allres_table)))
  colnames(freq_aa) <-  c("AAType", "VarType", "Frequency")
  freq_aa = subset(freq_aa, VarType %in% c("Cancer", "All_Residues") & !(AAType %in% c("UNK","SEC")), drop = T)
  freq_aa_sort <- freq_aa[with(freq_aa, order(Frequency)),]
  print(freq_aa_sort)
  freq_aa$AAType <- factor(freq_aa$AAType, levels = as.character(subset(freq_aa_sort, VarType == "Cancer")$AAType))
  print(freq_aa)
  
  p = ggplot(freq_aa, aes(x = AAType, fill = VarType)) + 
    geom_bar(aes(y = Frequency),stat="identity", position="dodge") + 
    scale_y_continuous(labels  = percent) + xlab("") + ylab("Frequency")
  ggsave(filename = "freq_aa_cosmic.pdf", height=3, width=12) 
}
freq_plot(protein_annotate_cancer_t)

freq_plot(protein_annotate_cancer_t)

######################## Allosteric on Cancer #####################################
allo_mix_cancer = sqldf("select allo_mix.*, protein_annotate_withcancer_u.VarType as cancer_var from allo_mix left join protein_annotate_withcancer_u")
fish_bs <- function(p_annotate_bs, vartype, loc){
  fish_result = fisher.test(table(data.frame(Disease = p_annotate_bs$cancer_var == vartype, Allosite = p_annotate_bs$allosite == loc)))
  fish_result
}
fish_bs(subset(allo_mix_cancer, allosite %in% c("Binding Site", "Allo")), "Cancer", "Allo")




############### aa change odds ratio #########
# turn out to be not right.
aa_change_preference <- function(p_annotate_bs, interested_par){
  aa_change <- factor(unique(p_annotate_bs[,interested_par]))
  vartypes   = c("Cancer")
  ratio_amino <- data.frame(merge(aa_change, vartypes, all = TRUE))
  colnames(ratio_amino) <- c(interested_par, "vartypes")
  for (vartype in vartypes) {
    for (each in aa_change){
      each = as.character(each)
      table_res <- table(deal_with_na(p_annotate_bs$VarType == vartype),
                         deal_with_na(as.character(p_annotate_bs[,interested_par]) == each))
      print(c(vartype, each))
      print(table_res)
      table_res <- apply(table_res, 1:2, as.numeric)
      rownames(table_res) <- c(paste("Not", vartype), vartype)
      colnames(table_res) <- c(paste("Not", each), each)
      
      print(table_res)
      
      fish_result = fisher.test(table(data.frame(residueName = (p_annotate_bs[,interested_par] == each), Disease = p_annotate_bs$VarType == vartype)))
      
      if (fish_result$estimate == 0) next
      
      ratio_amino[ratio_amino[,interested_par] == each & ratio_amino$vartypes == vartype,"estimate"] = fish_result$estimate
      ratio_amino[ratio_amino[,interested_par] == each & ratio_amino$vartypes == vartype,"ci_low"] = fish_result$conf.int[1]
      ratio_amino[ratio_amino[,interested_par] == each & ratio_amino$vartypes == vartype,"ci_up"] = fish_result$conf.int[2]
      ratio_amino[ratio_amino[,interested_par] == each & ratio_amino$vartypes == vartype,"pvalue"] = fish_result$p.value
    }
  }
  ratio_amino
}

protein_annotate_cancer_t[is.na(protein_annotate_cancer_t$VarType), "VarType"] = "Not Cancer"

ratio_amino <- aa_change_preference(protein_annotate_cancer_t, "uniprot_resnam_3d.before")

# for only protein surface, 
ratio_amino <- aa_change_preference(subset(protein_annotate_cancer_t, location == "Binding Site"), "uniprot_resnam_3d.before")
ratio_amino <- aa_change_preference(subset(protein_annotate_cancer_t, location == "Surface"), "uniprot_resnam_3d.before")


############################Plot odds ratio######################################
plot_ratio <- function(ratio_amino){
  prop = names(ratio_amino)[1]
  valid_var = unique(ratio_amino$vartypes)[1]
  ratio_amino_sort <- ratio_amino[with(ratio_amino, order(-estimate)),]
  ratio_amino[, prop] <- factor(ratio_amino[, prop], levels = as.character(subset(ratio_amino_sort, vartypes == valid_var)[, prop]))
  
  print(ratio_amino)
  ggplot(subset(ratio_amino, vartypes %in% c("Cancer", "Not Cancer")), aes_string(x = prop, y = 'estimate', ymin = 'ci_low', ymax = 'ci_up')) + geom_pointrange(aes(col = vartypes), position=position_dodge(width=0.30))  + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("") + scale_y_log10() + theme(axis.text.x = element_text( size=8, angle=30))
  ggsave(filename = "tmp.pdf", height=4, width=12)  
}
plot_ratio(subset(ratio_amino, !(uniprot_resnam_3d.before %in% c("SEC"))))



########################AA change after mutation##############################
idx <- with(protein_annotate_cancer_t, !is.na(uniprot_resnam_3d.before) & !is.na(uniprot_resnam_3d.after))
protein_annotate_cancer_change = protein_annotate_cancer_t
protein_annotate_cancer_change$aa_change <- NULL
protein_annotate_cancer_change[idx, "aa_change"] <- paste(protein_annotate_cancer_t[idx, 'uniprot_resnam_3d.before'], "to", protein_annotate_cancer_change[idx, 'uniprot_resnam_3d.after'])
ratio_amino <- aa_change_preference(protein_annotate_cancer_change, "aa_change")


############################Plot odds ratio######################################
plot_ratio <- function(ratio_amino){
  ratio_amino_sort <- ratio_amino[with(ratio_amino, order(-estimate)),]
  ratio_amino$aa_change <- factor(ratio_amino$aa_change, levels = as.character(subset(ratio_amino_sort, vartypes == "Cancer")$aa_change))
  ratio_amino <- subset(ratio_amino, pvalue < 0.05 & estimate > 0 & is.finite(estimate) & vartypes %in% c("Cancer", "Polymorphism"))
  common_change <- intersect(subset(ratio_amino, vartypes == "Cancer")$aa_change, subset(ratio_amino, vartypes == "Cancer")$aa_change)
  ratio_amino <- subset(ratio_amino, aa_change %in% common_change)
  
  ggplot(ratio_amino, aes(x = aa_change, y = estimate, ymin = ci_low, ymax = ci_up)) + geom_pointrange(aes(col = vartypes), position=position_dodge(width=0.30))  + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("") + scale_y_log10() + theme(axis.text.x = element_text(size=8, angle=30))
  ggsave(filename = "tmp.pdf", height=4, width=12)
}
plot_ratio(ratio_amino)