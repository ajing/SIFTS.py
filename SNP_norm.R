# Script for normalized SNP data set

protein_annotate<- read.table("./Result/avedist2surface.txt_small", sep = "\t", header = F)

# 
colnames(protein_annotate)<-c("pdbid", "chainid", "resnum", "resnam", "ssa", "rsa")

protein_annotate$location = "Core"

#MAX_ACC["GLY"] = 84.0
protein_annotate$location[protein_annotate$ssa > 5] = "Surface"

protein_annotate[protein_annotate$distance < 4.0 | !is.null(protein_annotate$BindingSiteComment),'location'] = "Binding Site"

table(protein_annotate$vartype, protein_annotate$location)

# The first table
allres_table <- table(protein_annotate$location)

get_stat_eachtype <- function(snp_type){
  loc_table = table(subset(protein_annotate, vartype == snp_type)$location)
  print("Observed")
  print(loc_table)
  print("Expected:")
  exp_table = allres_table / sum(allres_table) * sum(loc_table)
  print(exp_table)
  print("O/E")
  print(loc_table/exp_table)
}

library("abd")
library("fmsb")
oddsratio(x, conf.level = 0.95)

odds.ratio(x, conf.level = 0.95)
