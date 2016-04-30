# how about cancer mutation? does it also likely to be in protein ligand binding sites?

#################################################
#### For cosmic
#################################################
cosmic <- read.table("/users/ajing/SNPDist/Data/CosmicCompleteExport.tsv", sep = "\t", header = T, quote = "", na.string = "")
load("RData/cancer_u.RData")


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
protein_annotate_cancer[!is.na(protein_annotate_cancer$Mutation.AA.Before) & protein_annotate_cancer$FATHMM.prediction = "CANCER", "VarType"] = "Cancer"

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



allo_mix_cancer = sqldf("select allo_annotate.*, protein_annotate_withcancer_u.VarType as cancer_var from allo_annotate  left join protein_annotate_withcancer_u on allo_annotate.UniProtID = protein_annotate_withcancer_u.UniProtID and allo_annotate.uniprot_resnum = protein_annotate_withcancer_u.uniprot_resnum")
allo_mix_cancer_withsnp = subset(allo_mix_cancer, UniProtID %in% factor(unique(subset(allo_mix_cancer, !is.na(cancer_var))$UniProtID)))
allo_mix_cancer_withsnp$allosite = allo_mix_cancer_withsnp$location
allo_mix_cancer_withsnp$allosite[!is.na(allo_mix_cancer_withsnp$PubMedID)] = "Allo"

allo_mix_cancer_withsnp[is.na(allo_mix_cancer_withsnp$cancer_var),"cancer_var"] = "Not Cancer"


# binding site disease
unique(subset(protein_annotate_withsnp, location == "Binding Site")[,"DiseaseName"])  #182
# protein core disease
unique(subset(protein_annotate_withsnp, location == "Core")[,"DiseaseName"])  # 208

intersect(unique(subset(protein_annotate_withsnp, location == "Binding Site")[,"DiseaseName"]),unique(subset(protein_annotate_withsnp, location == "Core")[,"DiseaseName"]))   #115




#########
## 2/26/2016 change the definision of cancer SNPs, use FATHMM CANCER annotation
#########
library(data.table)
protein_annotate_withcancer_t <- data.table(protein_annotate_withcancer_u)
protein_annotate_withcancer_t[, VarType := "Passenger"] 
protein_annotate_withcancer_t[FATHMM.prediction == "CANCER", VarType := "Cancer"]

protein_annotate_withcancer_t[, table(VarType)]
get_stat_eachtype(protein_annotate_withcancer_t, "Cancer")
get_stat_eachtype(protein_annotate_withcancer_t, "Passenger")



# test the result by removing proteins with large number of mutations.
protein_annotate_withcancer_t <- data.table(protein_annotate_withcancer_u)
protein_annotate_withcancer_t[!is.na(Mutation.AA.Before), VarType := "Cancer"] 
num_snp_cancer <- protein_annotate_withcancer_t[!is.na(VarType), .N, by = UniProtID]

odds_table(protein_annotate_withcancer_t[!(UniProtID %in% num_snp_cancer[N>45, UniProtID])], "Cancer")
odds_table(protein_annotate_withcancer_t[!(UniProtID %in% num_snp_cancer[N>50, UniProtID])], "Cancer")
odds_table(protein_annotate_withcancer_t[!(UniProtID %in% num_snp_cancer[N>55, UniProtID])], "Cancer")
odds_table(protein_annotate_withcancer_t[!(UniProtID %in% num_snp_cancer[N>60, UniProtID])], "Cancer")
odds_table(protein_annotate_withcancer_t[!(UniProtID %in% num_snp_cancer[N>65, UniProtID])], "Cancer")
odds_table(protein_annotate_withcancer_t[!(UniProtID %in% num_snp_cancer[N>70, UniProtID])], "Cancer")


odds_table(protein_annotate_withcancer_t, "Cancer")

