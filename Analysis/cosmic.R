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
protein_annotate_cancer <- sqldf("select * from protein_annotate left join cosmic_missense where Gene.name = gene_name_acc and Mutation.AA.resnum = resnum")
protein_annotate_cancer$VarType = NA
protein_annotate_cancer[!is.na(protein_annotate_cancer$Mutation.AA.Before), "VarType"] = "Cancer"

table(protein_annotate_cancer$VarType, protein_annotate_cancer$location)

# The first table
get_stat_eachtype <- function(p_annotate, snp_type){
  allres_table <- table(p_annotate$location)
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
get_stat_eachtype(protein_annotate_cancer, "Cancer")


odds_ratio_stat <- function(p_annotate, vartype){
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
odds_ratio_stat(protein_annotate_withsnp, "Cancer")

# odds ratio between hydrophobic property
aa_prop_preference <- function(p_annotate_bs){
  prop_change <- factor(unique(p_annotate_bs$h_prop_change))
  vartypes   = c("Cancer")
  ratio_amino <- data.frame(merge(prop_change, vartypes, all = TRUE))
  colnames(ratio_amino) <- c("h_prop_change", "vartypes")
  for (vartype in vartypes) {
    for (each in prop_change){
      print(each)
      print(vartype)
      each = as.character(each)
      table_res <- table(deal_with_na(p_annotate_bs$VarType == vartype),
                         deal_with_na(as.character(p_annotate_bs$h_prop_change) == each))
      table_res <- apply(table_res, 1:2, as.numeric)
      rownames(table_res) <- c(paste("Not", vartype), vartype)
      colnames(table_res) <- c(paste("Not", each), each)
      
      fish_result = fisher.test(table(data.frame(residueName = (p_annotate_bs$h_prop_change == each), Disease = p_annotate_bs$VarType == vartype)))
      print(fish_result)
      print(ratio_amino$h_prop_change == each & ratio_amino$vartypes == vartype)
      ratio_amino[ratio_amino$h_prop_change == each & ratio_amino$vartypes == vartype,"estimate"] = fish_result$estimate
      ratio_amino[ratio_amino$h_prop_change == each & ratio_amino$vartypes == vartype,"ci_low"] = fish_result$conf.int[1]
      ratio_amino[ratio_amino$h_prop_change == each & ratio_amino$vartypes == vartype,"ci_up"] = fish_result$conf.int[2]
      ratio_amino[ratio_amino$h_prop_change == each & ratio_amino$vartypes == vartype,"pvalue"] = fish_result$p.value
    }
  }
  ratio_amino
}
# why protein_annotate_onlysnp, odds ratio for only changes
ratio_amino <- aa_prop_preference(protein_annotate_cancer)
