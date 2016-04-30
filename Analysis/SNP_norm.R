# Script for normalized SNP data set

protein_annotate<- read.table("./Data/SNPOnStruct_final.txt", sep = "\t", header = F, quote = "", na.string = "\\N")
# with aabefore and aaafter
protein_annotate<- read.table("./Data/SNPOnStruct_final_3.txt", sep = "\t", header = F, quote = "", na.string = "\\N")
# remove redundancy in gene level.
protein_annotate<- read.table("./Data/SNPOnStruct_final_7.txt", sep = "\t", header = F, quote = "", na.string = "\\N")
# add three letter amino acid type
protein_annotate<- read.table("./Data/SNPOnStruct_final_8.txt", sep = "\t", header = F, quote = "", na.string = "\\N")


#
colnames(protein_annotate)<-c("pdbid", "biounit", "ModelID", "chainid", "resnam", "resnum", "structcode", "ssa", "rsa", "UniProtID", "uniprot_resnam", "uniprot_resnum", "bs_biounit", "bs_p_chainid", "bs_p_resnum", "ligandName", "BindingSiteComment", "distance", "genename", "SwissProt_AC", "FTID", "AABefore", "prot_resnum", "AAAfter", "VarType", "dbSNPID", "DiseaseName", "uniprot_resnam_3d", "AABeforeProp", "AAAfterProp", "proteinname", "reviewed", "gene_name_acc", "interproname_acc")

protein_annotate$location = "Core"

#MAX_ACC["GLY"] = 84.0
protein_annotate$location[protein_annotate$ssa > 5] = "Surface"

protein_annotate[ which(protein_annotate$distance <= 4.0 | !is.na(protein_annotate$BindingSiteComment)),'location'] = "Binding Site"

table(protein_annotate$VarType, protein_annotate$location)


# only keep UniProt with at least one nsSNP
snp_num = ddply(protein_annotate, .(UniProtID), summarise, N = length(unique(uniprot_resnum[!is.na(VarType)])))
protein_annotate_withsnp = subset(protein_annotate, UniProtID %in% factor(unique(subset(protein_annotate, !is.na(VarType))$UniProtID)))
length(unique(protein_annotate_withsnp$UniProtID))

> table(protein_annotate_withsnp$location)

Binding Site         Core      Surface 
17010        49292       154993 

> length(unique(protein_annotate[protein_annotate$VarType == "Unclassified",]$UniProtID))
[1] 150
> length(unique(protein_annotate[protein_annotate$VarType == "Disease",]$UniProtID))
[1] 246
> length(unique(protein_annotate[protein_annotate$VarType == "Polymorphism",]$UniProtID))
[1] 581


> dim(unique(protein_annotate[protein_annotate$VarType == "Polymorphism",][,c("UniProtID", "uniprot_r
esnum")]))
[1] 1656    2
> dim(unique(protein_annotate[protein_annotate$VarType == "Disease",][,c("UniProtID", "uniprot_resnum
")]))
[1] 3006    2
> dim(unique(protein_annotate[protein_annotate$VarType == "Unclassified",][,c("UniProtID", "uniprot_r
esnum")]))
[1] 813   2


> dim(unique(protein_annotate[protein_annotate$VarType == "Unclassified",]))
[1] 1766   34
> dim(unique(protein_annotate[protein_annotate$VarType == "Polymorphism",]))
[1] 1717   34
> dim(unique(protein_annotate[protein_annotate$VarType == "Disease",]))
[1] 3709   34



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
get_stat_eachtype(protein_annotate_withsnp, "Disease")

ggplot(protein_annotate, aes(x=ssa)) +
  geom_histogram() + facet_grid(.~location)

ggplot(subset(protein_annotate, location %in% c("Surface", "Binding Site")), aes(x=ssa, color = location)) +
  geom_density() + facet_grid(location~.)

ggsave(filename = "tmp.pdf")

summary(subset(protein_annotate, location == "Binding Site" & VarType =="Polymorphism"))

# odds ratio is to quantify how strongly the presence or absence of property A is associated with the presence or absense of property B in a given population.
# If the OR is greater than 1, then having "A" is considered to be "associated" with having "B" in the sense that the having of "B" raises (relative to not-having "B") the odds of having "A".

library("fmsb")
result = oddsratio(, conf.level = 0.95)
result$p.value
result$estimate
result$conf.int

install.packages(epitools)
library("epitools")

deal_with_na <- function(x) {
  x[is.na(x)] = F
  factor(x)
}

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

# Disease
odds_ratio_stat(protein_annotate_withsnp, "Disease")

# Polymorphism
odds_ratio_stat(protein_annotate_withsnp, "Polymorphism")

# Polymorphism without antigen
odds_ratio_stat(subset(protein_annotate_withsnp, !grepl("*histocompatibility antigen*", x = proteinname)), "Polymorphism")
with(subset(protein_annotate_withsnp, !grepl("*histocompatibility antigen*", x = proteinname)), table(location, vartype))


# Unclassified
odds_ratio_stat(protein_annotate_withsnp, "Unclassified")

# Fisher's Exact Test
# examine the significance of the association (contingency) between the two kinds of classification.
# null hypothesis that disease and non-disease are equally likely to be in the structure core,
fish_result = fisher.test(table(DF$Admit, DF$Gender))
fish_result$p.value
fish_result$estimate
fish_result$conf.int

# For amino acid mutated in binding site
# the subset I am interested in is subset(protein_annotate, 'location' == "Binding Site")
p_annotate_bs <- subset(protein_annotate, location == "Binding Site")
p_annotate_bs <- subset(protein_annotate, location %in% c("Surface", "Binding Site"))
aa_preference <- function(p_annotate_bs){
  amino_names <- factor(unique(p_annotate_bs$resnam)[1:20])
  each = as.character(amino_names[1])
  vartype = "Disease"
  vartypes   = c("Disease", "Polymorphism", "Unclassified")
  ratio_amino <- data.frame(merge(amino_names, vartypes, all = TRUE))
  colnames(ratio_amino) <- c("amino_names", "vartypes")
  for (vartype in vartypes) {
    for (each in amino_names){
      #print(each)
      #print(vartype)
      each = as.character(each)
      table_res <- table(deal_with_na(p_annotate_bs$VarType == vartype),
                         deal_with_na(as.character(p_annotate_bs$resnam) == each))
      table_res <- apply(table_res, 1:2, as.numeric)
      rownames(table_res) <- c(paste("Not", vartype), vartype)
      colnames(table_res) <- c(paste("Not", each), each)
      
      fish_result = fisher.test(table(data.frame(residueName = (p_annotate_bs$resnam == each), Disease = p_annotate_bs$VarType == vartype)))
      print(fish_result)
      print(ratio_amino$amino_names == each & ratio_amino$vartypes == vartype)
      ratio_amino[ratio_amino$amino_names == each & ratio_amino$vartypes == vartype,"estimate"] = fish_result$estimate
      ratio_amino[ratio_amino$amino_names == each & ratio_amino$vartypes == vartype,"ci_low"] = fish_result$conf.int[1]
      ratio_amino[ratio_amino$amino_names == each & ratio_amino$vartypes == vartype,"ci_up"] = fish_result$conf.int[2]
      ratio_amino[ratio_amino$amino_names == each & ratio_amino$vartypes == vartype,"pvalue"] = fish_result$p.value
    }
  }
  ratio_amino
}
ratio_amino <- aa_preference(protein_annotate_withsnp)
ratio_amino <- aa_preference(subset(protein_annotate_withsnp, location %in% c("Surface", "Binding Site")))
ratio_amino <- aa_preference(subset(protein_annotate_withsnp, location %in% c("Core")))
ratio_amino <- aa_preference(subset(protein_annotate_withsnp, location %in% c("Binding Site")))

# write to table
write.csv(ratio_amino, file = "tmp.csv")

# for the amino acid it mutated to
protein_annotate_withsnp_cp <- protein_annotate_withsnp
protein_annotate_withsnp_cp$resnam[!is.na(protein_annotate_withsnp$AAAfter)] <- protein_annotate_withsnp$AAAfter[!is.na(protein_annotate_withsnp$AAAfter)]

ratio_amino <- aa_preference(subset(protein_annotate_withsnp_cp, location %in% c("Surface", "Binding Site")))

############################Plot odds ratio######################################
ratio_amino_sort <- ratio_amino[with(ratio_amino, order(-mean)),]
ratio_amino$amino_names <- factor(ratio_amino$amino_names, levels = as.character(subset(ratio_amino_sort, vartypes == "Disease")$amino_names))

ggplot(subset(ratio_amino, vartypes %in% c("Disease", "Polymorphism")), aes(x = amino_names, y = estimate, ymin = ci_low, ymax = ci_up)) + geom_pointrange(aes(col = vartypes), position=position_dodge(width=0.30))  + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("") + scale_y_log10()
ggsave(filename = "tmp.pdf", height=3, width=12) 



# odds ratio between binding sites and other surface residues
aa_loc_preference <- function(p_annotate_bs){
  amino_names <- factor(unique(p_annotate_bs$resnam)[1:20])
  each = as.character(amino_names[1])
  locs   = c("Binding Site")
  ratio_amino <- data.frame(merge(amino_names, locs, all = TRUE))
  colnames(ratio_amino) <- c("amino_names", "location")
  for (loc in locs) {
    for (each in amino_names){
      #print(each)
      each = as.character(each)
      table_res <- table(deal_with_na(p_annotate_bs$location == loc),
                         deal_with_na(as.character(p_annotate_bs$resnam) == each))
      table_res <- apply(table_res, 1:2, as.numeric)
      print(table_res)
      rownames(table_res) <- c(paste("Not", loc), loc)
      colnames(table_res) <- c(paste("Not", each), each)
      
      fish_result = fisher.test(table(data.frame(residueName = (p_annotate_bs$resnam == each), Disease = p_annotate_bs$location == loc)))
      print(fish_result)
      print(ratio_amino$amino_names == each & ratio_amino$location == loc)
      ratio_amino[ratio_amino$amino_names == each & ratio_amino$location == loc,"estimate"] = fish_result$estimate
      ratio_amino[ratio_amino$amino_names == each & ratio_amino$location == loc,"ci_low"] = fish_result$conf.int[1]
      ratio_amino[ratio_amino$amino_names == each & ratio_amino$location == loc,"ci_up"] = fish_result$conf.int[2]
      ratio_amino[ratio_amino$amino_names == each & ratio_amino$location == loc,"pvalue"] = fish_result$p.value
    }
  }
  ratio_amino
}
ratio_amino <- aa_loc_preference(subset(protein_annotate_withsnp, location %in% c("Surface", "Binding Site")))

ratio_amino_sort <- ratio_amino[with(ratio_amino, order(-estimate)),]
ratio_amino$amino_names <- factor(ratio_amino$amino_names, levels = as.character(subset(ratio_amino_sort, location == "Binding Site")$amino_names))
ggplot(ratio_amino, aes(x = amino_names, y = estimate, ymin = ci_low, ymax = ci_up)) + geom_pointrange(color = I("blue"))  + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("") + scale_y_log10()
ggsave(filename = "tmp.pdf", height=3, width=12) 


# how many disease covered
> length(unique(protein_annotate_withsnp$DiseaseName))
[1] 380


# pearson correlation coefficient
cor(subset(ratio_amino, vartypes == "Disease")$mean, subset(ratio_amino, vartypes == "Polymorphism")$mean, method="pearson") 
cor.test(subset(ratio_amino, vartypes == "Disease")$mean, subset(ratio_amino, vartypes == "Polymorphism")$mean, method="pearson") 


# plots odds ratio with error bar
ggplot(dfwc.between, aes(x=condition, y=value, group=1)) +
  geom_bar() +
  geom_errorbar(width=.1, aes(ymin=value-ci, ymax=value+ci), colour="red")


library(ggplot2)
library(reshape2)
library(scales)
# frequency distribution
freq_aa = melt(cbind(prop.table(table(protein_annotate_withsnp$uniprot_resnam_3d, protein_annotate_withsnp$VarType), margin = 2),All_Residues = prop.table(table(protein_annotate_withsnp$uniprot_resnam_3d))))
colnames(freq_aa) <-  c("AAType", "VarType", "Frequency")
freq_aa = subset(freq_aa, VarType %in% c("Disease","Polymorphism", "All_Residues") & !(AAType %in% c("UNK","SEC")), drop = T)
freq_aa_sort <- freq_aa[with(freq_aa, order(Frequency)),]
freq_aa$AAType <- factor(freq_aa$AAType, levels = as.character(subset(freq_aa_sort, VarType == "Disease")$AAType))

p = ggplot(freq_aa, aes(x = AAType, fill = VarType)) + 
  geom_bar(aes(y = Frequency),stat="identity", position="dodge") + 
  scale_y_continuous(labels  = percent) + xlab("") + ylab("Frequency")
ggsave(filename = "freq_aa.pdf", height=3, width=12) 

# freq plot for location
freq_loc <- melt(prop.table(table(protein_annotate_withsnp$uniprot_resnam_3d, protein_annotate_withsnp$location), margin = 2))
colnames(freq_loc) <-  c("AAType", "Location", "Frequency")
freq_loc = subset(freq_loc, !(AAType %in% c("UNK","SEC","NA")), drop = T)
freq_loc_sort <- freq_aa[with(freq_loc, order(Frequency)),]
freq_loc$AAType <- factor(freq_loc$AAType, levels = as.character(subset(freq_loc_sort, VarType == "Disease")$AAType))

p = ggplot(freq_loc, aes(x = AAType, fill = Location)) + 
  geom_bar(aes(y = Frequency),stat="identity", position="dodge") + 
  scale_y_continuous(labels  = percent) + xlab("") + ylab("Frequency")
ggsave(filename = "freq_loc.pdf", height=3, width=12) 


p = ggplot(subset(protein_annotate, VarType == "Disease" & !(uniprot_resnam_3d %in% c("UNK","SEC")), drop = T), aes(x = uniprot_resnam_3d)) + 
  geom_bar(aes(y = (..count..))) + xlab("") + ylab("Frequency")
ggsave(filename = "freq_aa.pdf", height=9, width=12) 

# http://stackoverflow.com/questions/13386177/how-to-create-odds-ratio-and-95-ci-plot-in-r
# http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_%28ggplot2%29/
# http://stackoverflow.com/questions/14069629/plotting-confidence-intervals
ggplot(subset(ratio_amino, vartypes %in% c("Disease", "Polymorphism")), aes(x = amino_names, y = mean, ymin = ci_low, ymax = ci_up)) + geom_pointrange(aes(col = vartypes), position=position_dodge(width=0.30))  + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("") + scale_y_log10()
ggsave(filename = "tmp.pdf", height=9, width=12) 

# Also odds ratio for each type of amino acids
# for disease and polymorphism seperately
ggplot(dat, aes(x = pollut, y = or, ymin = lcl, ymax = ucl)) + geom_pointrange(aes(col = factor(lag)), position=position_dodge(width=0.30)) + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("")

# odds ratio comparison between disease and polymorphism
ggplot(dat, aes(x = pollut, y = or, ymin = lcl, ymax = ucl)) + geom_pointrange(aes(col = factor(lag)), position=position_dodge(width=0.30)) + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("")


# Distribution of BLOSSUM Matrix
# 3x3 chis-sq test







# Property change after mutation
protein_annotate_withsnp$type = paste(protein_annotate_withsnp$AABeforeProp, "to", protein_annotate_withsnp$AAAfterProp)
protein_annotate_onlysnp = subset(protein_annotate_withsnp, !is.na(VarType))
# PH = 2
hydro_prop <- data.frame(AAName = c("LEU", "ILE", "PHE", "TRP", "VAL", "MET", "CYS", "TYR", "ALA", "THR", "GLU", "GLY", "SER", "GLN", "ASP", "ARG", "LYS", "ASN", "HIS", "PRO"), H_Prop = c("Very Hydrophobic", "Very Hydrophobic", "Very Hydrophobic", "Very Hydrophobic", "Very Hydrophobic", "Very Hydrophobic", "Hydrophobic", "Hydrophobic", "Hydrophobic", "Neutral", "Neutral","Neutral","Neutral","Neutral","Neutral", "Hydrophilic", "Hydrophilic", "Hydrophilic","Hydrophilic","Hydrophilic"))
# PH = 7
hydro_prop <- data.frame(AAName = c("LEU", "ILE", "PHE", "TRP", "VAL", "MET", "CYS", "TYR", "ALA", "THR", "GLU", "GLY", "SER", "GLN", "ASP", "ARG", "LYS", "ASN", "HIS", "PRO"), H_Prop = c("Very Hydrophobic", "Very Hydrophobic", "Very Hydrophobic", "Very Hydrophobic", "Very Hydrophobic", "Very Hydrophobic", "Hydrophobic", "Hydrophobic", "Hydrophobic", "Neutral", "Neutral","Neutral","Neutral","Neutral","Hydrophilic", "Hydrophilic", "Hydrophilic", "Hydrophilic","Hydrophilic","Hydrophilic"))

#protein_annotate_withsnp$h_prop_change <- apply(protein_annotate_withsnp, 1, function(x) { paste(as.character(hydro_prop[hydro_prop$AAName == x[["AABefore"]], "H_Prop"]), "to", as.character(hydro_prop[hydro_prop$AAName == x[["AAAfter"]], "H_Prop"]))})

protein_annotate_onlysnp$h_prop_change <- apply(protein_annotate_onlysnp, 1, function(x) { paste(as.character(hydro_prop[hydro_prop$AAName == x[["AABefore"]], "H_Prop"]), "to", as.character(hydro_prop[hydro_prop$AAName == x[["AAAfter"]], "H_Prop"]))})


# odds ratio between hydrophobic property
aa_prop_preference <- function(p_annotate_bs){
  prop_change <- factor(unique(p_annotate_bs$h_prop_change))
  vartype = "Disease"
  vartypes   = c("Disease", "Polymorphism", "Unclassified")
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
ratio_amino <- aa_prop_preference(protein_annotate_onlysnp)


############################Plot odds ratio######################################
ratio_amino_sort <- ratio_amino[with(ratio_amino, order(-estimate)),]
ratio_amino$h_prop_change <- factor(ratio_amino$h_prop_change, levels = as.character(subset(ratio_amino_sort, vartypes == "Disease")$h_prop_change))

ggplot(subset(ratio_amino, vartypes %in% c("Disease", "Polymorphism")), aes(x = h_prop_change, y = estimate, ymin = ci_low, ymax = ci_up)) + geom_pointrange(aes(col = vartypes), position=position_dodge(width=0.30))  + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("") + scale_y_log10() + theme(axis.text.x = element_text( size=8, angle=30))
ggsave(filename = "tmp.pdf", height=4, width=12) 


# continue analysis the hydrophibicity
with(subset(protein_annotate_onlysnp,h_prop_change == "Hydrophobic to Hydrophobic"), table(factor(VarType), factor(AABefore), factor(AAAfter)))

# For hydrophobic to hydrophobic ranking first because all of them are related to mutation from Cysteine to Tyrosine (72 cases) and Tyrosine to Cysteine (60 cases). From CYS to TYR, the protein will lose disulfide bond. From TYR to CYS, the protein will lose an aromatic ring. (add the codon mutate from CYS to TYR and TYR to CYS)


########################AA change after mutation##############################
protein_annotate_onlysnp$aa_change <- paste(protein_annotate_onlysnp$AABefore, "to", protein_annotate_onlysnp$AAAfter)
# odds ratio between amino acid change
aa_change_preference <- function(p_annotate_bs, interested_par){
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
ratio_amino <- aa_change_preference(protein_annotate_onlysnp, "aa_change")

############################Plot odds ratio######################################
ratio_amino_sort <- ratio_amino[with(ratio_amino, order(-estimate)),]
ratio_amino$aa_change <- factor(ratio_amino$aa_change, levels = as.character(subset(ratio_amino_sort, vartypes == "Disease")$aa_change))
ratio_amino <- subset(ratio_amino, pvalue < 0.05 & estimate > 0 & is.finite(estimate) & vartypes %in% c("Disease", "Polymorphism"))
common_change <- intersect(subset(ratio_amino, vartypes == "Disease")$aa_change, subset(ratio_amino, vartypes == "Polymorphism")$aa_change)
ratio_amino <- subset(ratio_amino, aa_change %in% common_change)

ggplot(ratio_amino, aes(x = aa_change, y = estimate, ymin = ci_low, ymax = ci_up)) + geom_pointrange(aes(col = vartypes), position=position_dodge(width=0.30))  + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("") + scale_y_log10() + theme(axis.text.x = element_text(size=8, angle=30))
ggsave(filename = "tmp.pdf", height=4, width=12) 






########## Boostrap


########## TYR to CYS or CYS to TYR
# side chain size size change CYS(135) to TYR(222) and disulphide bond
> with(subset(protein_annotate_onlysnp, aa_change=="TYR to CYS"), length(table(factor(DiseaseName)))) 
[1] 56
> with(subset(protein_annotate_onlysnp, aa_change=="CYS to TYR"), length(table(factor(DiseaseName)))) 
[1] 44

> with(subset(protein_annotate_onlysnp, aa_change=="TYR to CYS"), length(table(factor(proteinname)))) 
[1] 58
> with(subset(protein_annotate_onlysnp, aa_change=="CYS to TYR"), length(table(factor(proteinname)))) 
[1] 44


########## molecular weight / size of side chain
sc_v = c(135.0, 163.0, 130.0, 157.0, 198.0, 205.0, 169.0, 136.0, 142.0, 197.0, 106.0, 84.0, 184.0, 164.0, 248.0, 227.0, 142.0, 194.0, 222.0, 188.0)
#hist(sc_v)
#sum(sc_v < 170)
sc_volumn = data.frame(AAName = c("ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"), volumn = c(106.0,135.0,163.0,194.0,197.0,84.0,184.0,169.0,205.0,164.0,188.0,157.0,136.0,198.0,248.0,130.0,142.0,142.0,227.0,222.0))
sc_volumn$volumn.c = cut(sc_volumn$volumn, breaks=3, label = c("small", "medium", "large"))
protein_annotate_onlysnp$volumn_change <- apply(protein_annotate_onlysnp, 1, function(x) { paste(as.character(sc_volumn[sc_volumn$AAName == x[["AABefore"]], "volumn.c"]), "to", as.character(sc_volumn[sc_volumn$AAName == x[["AAAfter"]], "volumn.c"]))})
ratio_amino <- aa_change_preference(protein_annotate_onlysnp, "volumn_change")

# save sc_volumn to file
write.csv(sc_volumn[order(sc_volumn$volumn),], file = "tmp.csv")

# for only snps in binding site
ratio_amino <- aa_change_preference(subset(protein_annotate_onlysnp, location == "Binding Site"), "volumn_change")

############################Plot odds ratio######################################
ratio_amino_sort <- ratio_amino[with(ratio_amino, order(-estimate)),]
ratio_amino$volumn_change <- factor(ratio_amino$volumn_change, levels = as.character(subset(ratio_amino_sort, vartypes == "Disease")$volumn_change))

ggplot(subset(ratio_amino, vartypes %in% c("Disease", "Polymorphism")), aes(x = volumn_change, y = estimate, ymin = ci_low, ymax = ci_up)) + geom_pointrange(aes(col = vartypes), position=position_dodge(width=0.30))  + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("") + scale_y_log10() + theme(axis.text.x = element_text( size=10))
ggsave(filename = "tmp.pdf", height=4, width=12) 


########### more analysis about size
> with(subset(protein_annotate_onlysnp, volumn_change %in% c("large to small", "small to large")), table(aa_change))
aa_change
ALA to GLU ALA to PHE ARG to CYS ARG to GLY ARG to PRO ARG to SER CYS to ARG 
25          3        149         87         73         42         56 
CYS to PHE CYS to TRP CYS to TYR GLN to PRO GLU to ALA GLU to GLY GLU to PRO 
30         19         75         19         32         50          1 
GLY to ARG GLY to GLU GLY to LYS GLY to PHE GLY to TRP LYS to GLY PHE to CYS 
186         82          1          1         14          1         28 
PHE to GLY PHE to SER PRO to ARG PRO to GLN PRO to PHE SER to ARG SER to LYS 
1         51         54          9          5         46          1 
SER to PHE SER to TRP SER to TYR TRP to CYS TRP to GLY TRP to SER TYR to CYS 
52          5         19         19         13          9         88 
TYR to GLY TYR to SER 
1         21 




################################## Blosum62 analysis##############################
library(seqinr)
library(peplib)
library(reshape2)

data(blosum62)
blosum_melt<-melt(blosum62)
blosum_melt <- subset(blosum_melt, !(Var1 %in% c("GAP", "B", "Z", "X")) & !(Var2 %in% c("GAP", "B", "Z", "X")))
blosum_melt$volumn_change <- apply(blosum_melt, 1, function(x) { alist = c(as.character(sc_volumn[sc_volumn$AAName == toupper(aaa(x[["Var1"]])), "volumn.c"]), as.character(sc_volumn[sc_volumn$AAName == toupper(aaa(x[["Var2"]])), "volumn.c"])); alist = sort(alist); print(alist); paste(alist[1], "to", alist[2])})
blosum_melt$value.mean = apply(blosum_melt, 1, function(x){mean(subset(blosum_melt, volumn_change==x[["volumn_change"]])$value)})

ggplot(blosum_melt, aes(x=value)) + geom_histogram() + geom_vline(aes(xintercept=value.mean),linetype="dashed", size=1, color = "red") + facet_grid(volumn_change~.)


########################## Get the table after removing the antigen protein
no_antigen_table = with(subset(protein_annotate_withsnp, !grepl("*histocompatibility antigen*", x = proteinname)), table(location, VarType))

########################## Linear Model
# the size of side chain sc_volumn
# the hydrophobic porperty, 
hydro_index <- data.frame(AAName = c("PHE", "ILE", "TRP", "LEU", "VAL", "MET", "TYR", "CYS", "ALA", "THR", "HIS", "GLY", "SER", "GLN", "ARG", "LYS", "ASN", "GLU", "PRO", "ASP"), index_num = c( 100, 99, 97, 97, 76, 74, 63, 49, 41, 13, 8, 0, -5, -10, -14, -23, -28, -31, -46, -55))
protein_annotate_onlysnp$hydro_change <- apply(protein_annotate_onlysnp, 1, function(x) { hydro_index[hydro_index$AAName == x[["AABefore"]], "index_num"] - hydro_index[hydro_index$AAName == x[["AAAfter"]], "index_num"]})
protein_annotate_onlysnp$size_change <- apply(protein_annotate_onlysnp, 1, function(x) { sc_volumn[sc_volumn$AAName == x[["AABefore"]], "volumn"] - sc_volumn[sc_volumn$AAName == x[["AAAfter"]], "volumn"]})
protein_annotate_onlysnp$is_disease <- protein_annotate_onlysnp$VarType == "Disease"

model <- lm(is_disease ~ hydro_change + location + location:size_change, data = protein_annotate_onlysnp)

#model <- lm(is_disease ~ hydro_change:loc_factor, data = protein_annotate_onlysnp)
summary(model)
pre_result_lm <- predict(model)


library('e1071')
model <- svm(is_disease ~ ., data = subset(protein_annotate_onlysnp, select = c("is_disease", "hydro_change", "size_change", "location")), type = "C-classification", probability = T)
pre_result <- predict(model, subset(protein_annotate_onlysnp, select = c( "hydro_change", "size_change", "location")), probability = T)
pre_result_svm <- attr(pre_result, "probabilities")[,2]

# load data from SIFT
SIFT_score <- read.table("./Data/siftscore.txt", sep = "\t", header = F, quote = "", na.string = "\\N")
colnames(SIFT_score) <- c("FTID", "SIFT_score")
SIFT_result<- subset(merge(protein_annotate_onlysnp, SIFT_score, by = c("FTID")), select = c("is_disease", "SIFT_score"))

pred <- with(SIFT_result, prediction(SIFT_score, is_disease))
perf_sift <- performance(pred, measure = "tpr", x.measure = "fpr")
auc_sift <- performance(pred, measure = "auc")
perf_sift$type <- "SIFT"
pdf("tmp.pdf")
plot(perf, col=rainbow(10))
abline(a = 0, b = 1)
dev.off()

#AUC
get_auc <- function(pre_result){
  pred <- prediction(pre_result, protein_annotate_onlysnp$is_disease)
  perf_loc <- performance(pred, measure = "tpr", x.measure = "fpr")
  auc_loc <- performance(pred, measure = "auc")
  attr(auc_loc, "y.value")[[1]]
}
get_auc(pre_result_lm)

install.packages("ROCR")
library(ROCR)

roc_d <- function(pre_result) {
  pred <- prediction(pre_result, protein_annotate_onlysnp$is_disease)
  perf <- performance(pred, measure = "tpr", x.measure = "fpr")
  data.frame(fpr = perf@x.values[[1]], tpr = perf@y.values[[1]])
}
roc_data <- rbind(cbind(roc_d(pre_result_lm), predictor = "Logistic"), cbind(roc_d(pre_result_svm), predictor = "SVM"), cbind(fpr = perf_sift@x.values[[1]], tpr = perf_sift@y.values[[1]], predictor = "SIFT"))

auc_loc <- performance(pred, measure = "auc")
perf_loc$type <- "location"
pdf("tmp.pdf")
plot(perf, col=rainbow(10))
abline(a = 0, b = 1)
dev.off()



########################## Gene enrichment analysis
get_gene_list <- function(genelist) {
  uniq_gene_1 <- unique(genelist)
  uniq_gene_2 <- c()
  for (i in 1:length(uniq_gene_1)) {
    uniq_gene_2 = c(uniq_gene_2, strsplit(as.character(uniq_gene_1[i]),',',fixed=TRUE)[[1]])
  }
  uniq_gene_2 = unique(as.character(uniq_gene_2)[!is.na(uniq_gene_2)])
}
gene_disease <- get_gene_list(subset(protein_annotate_withsnp, VarType == "Disease")$gene_name_acc)
gene_disease_bs <- get_gene_list(subset(protein_annotate_withsnp, VarType == "Disease" & location == "Binding Site")$gene_name_acc)
gene_disease_core <- get_gene_list(subset(protein_annotate_withsnp, VarType == "Disease" & location == "Core")$gene_name_acc)
gene_all <- get_gene_list(protein_annotate_withsnp$gene_name_acc)

gene_polymor <- get_gene_list(subset(protein_annotate_withsnp, VarType == "Polymorphism")$gene_name_acc)
gene_polymor_bs <- get_gene_list(subset(protein_annotate_withsnp, VarType == "Polymorphism" & location == "Binding Site")$gene_name_acc)

write(gene_all, file = "gene_all.txt", sep = "\n")
write(gene_disease, file = "gene_disease.txt", sep = "\n")
write(gene_disease_bs, file = "gene_disease_bs.txt", sep = "\n")
write(gene_disease_core, file = "gene_disease_core.txt", sep = "\n")

write(gene_polymor,  file = "gene_polymor.txt", sep = "\n")
write(gene_polymor_bs,  file = "gene_polymor_bs.txt", sep = "\n")

write(gene_polymor_bs,  file = "gene_polymor_bs.txt", sep = "\n")



################ more anlaysis for drug resistant

no_antigen_table = with(subset(protein_annotate_withsnp, !grepl("*transporter*", x = proteinname)), table(location, VarType))


################ for popular ligand in disease SNP binding site
tmp <- table(factor(subset(protein_annotate_onlysnp, location == "Binding Site" & VarType == "Disease", select = c("ligandName"))$ligandName))
sort(tmp)


################  Drug related information
drug_gene = read.table("./Data/drugbank_genelist.csv", header = F)
dim(subset(protein_annotate_onlysnp, location == "Binding Site" & VarType == "Disease" & genename %in% drug_gene))



################ Graph based analysis
ligand_smile = read.csv("./Data/ligand_smile.txt", sep="," , header = F)
names(ligand_smile) <- c("ligname", "smiles")

unique_ligand <- count(data.frame(subset(protein_annotate_onlysnp, !is.na(ligandName), select = c("VarType", "ligandName", "location"))))
unique_ligand <- merge(unique_ligand, ligand_smile, by.x = "ligandName", by.y = "ligname")
write.table(unique_ligand, file = "./Data/ligand_net.txt", sep="\t", quote = F, row.names = F)


protein_ligand <- count(data.frame(subset(protein_annotate_onlysnp, !is.na(ligandName), select = c("UniProtID", "ligandName"))))
subset(protein_annotate_onlysnp, VarType == "Disease" & !is.na(ligandName), select = c("UniProtID", "ligandName"))
write.table(protein_ligand, file = "./Data/prolig_net.txt")

poligand = read.csv("./Data/poly_ligand.txt", sep="\t" , header = F)
names(poligand) <- c("UniProtID", "uniprot_resnum", "VarType", "ligcount")
poligand$ligcount.mean = apply(poligand, 1, function(x){mean(subset(poligand, VarType==x[["VarType"]])$ligcount)})
poligand$total_res = apply(poligand, 1, function(x){length(subset(poligand, VarType==x[["VarType"]])$ligcount)})


ggplot(poligand, aes(x = ligcount, color = VarType)) + scale_x_log10() + geom_density() + geom_vline(aes(xintercept=ligcount.mean),linetype="dashed", size=1, color = "red") 
ggsave(filename = "tmp.pdf")


pdf("tmp.pdf")
hist(poligand$ligcount)
dev.off()

pdf("tmp.pdf")
with(poligand, boxplot(ligcount~VarType,notch=T, log = "y"))
dev.off()

ggplot(data=poligand, aes(x=VarType, y=ligcount, fill=VarType))  + scale_fill_discrete(guide = FALSE) + geom_boxplot(notch = T) + geom_point(aes(y = ligcount.mean), shape = 18, size = 3) + geom_text(aes(label = paste("mean = ", round(ligcount.mean, digits = 2), ", n = ", total_res, sep = ""), y = ligcount.mean + 1)) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + annotation_logticks(sides = "l") + xlab("SNV type") + ylab("Number of unique ligands")
ggsave(filename = "tmp.pdf")

with(poligand, wilcox.test(ligcount[VarType == "Disease"], ligcount[VarType == "Polymorphism"]))
with(poligand, wilcox.test(ligcount[VarType == "Polymorphism"], ligcount[VarType == "Unclassified"]))



####################### Central residues and peripheral residue
protein_annotate_withsnp_site <- subset(protein_annotate_withsnp, location == "Binding Site")
pro_lig_count <- read.table("./Data/pro_lig_count.txt", sep = "\t", header = F, quote = "", na.string = "\\N")
pro_res_lig <- read.table("./Data/pro_res_count_withsnp.txt", sep = "\t", header = F, quote = "", na.string = "\\N")
colnames(pro_lig_count) <- c("UniProtID", "ligcount")
colnames(pro_res_lig) <- c("UniProtID", "uniprot_resnum", "ligcount_res")
protein_annotate_withsnp_site <- merge(protein_annotate_withsnp_site, pro_lig_count, all.x = TRUE)
protein_annotate_withsnp_site <- merge(protein_annotate_withsnp_site, pro_res_lig, all.x = TRUE)

#protein_annotate_withsnp_site$site_annotat <- with(protein_annotate_withsnp_site, c("center", "periphery")[ 1 + (ligcount_res / ligcount < .75)])
protein_annotate_withsnp_site$site_annotat <- with(protein_annotate_withsnp_site, c("center", "periphery")[ 1 + (ligcount_res / ligcount < .9)])

fish_bs <- function(p_annotate_bs, vartype, loc){
  fish_result = fisher.test(table(data.frame(Disease = p_annotate_bs$VarType == vartype, Central = p_annotate_bs$site_annotat == loc)))
  fish_result
}
fish_bs(protein_annotate_withsnp_site, "Disease", "center")
fish_bs(protein_annotate_withsnp_site, "Polymorphism", "center")
fish_bs(protein_annotate_withsnp_site, "Unclassified", "center")
# based on this result polymorphism are more likely to be in the center of binding sites, and unclassifed SNP are less likely to be located in the center and disesae SNP has no preference

fish_bs(subset(protein_annotate_withsnp_site, !grepl("*histocompatibility antigen*", x = proteinname)), "Polymorphism", "center")




################# allsoteric sites or other sites
allo_res <- read.table("./Data/pro_res_count.txt", sep = "\t", header = F, quote = "", na.string = "\\N")
colnames(allo_res) <- c("UniProtID", "uniprot_resnum")
allo_res$allosite = "Allo"
allo_site <- subset(protein_annotate_withsnp, as.character(UniProtID) %in% as.character(unique(allo_res$UniProtID)))
allo_site <- merge(allo_site, allo_res, all.x = T)
allo_site[is.na(allo_site$allosite), "allosite"] <- allo_site[is.na(allo_site$allosite), "location"]
allo_site$allosite <- factor(allo_site$allosite)
fish_bs <- function(p_annotate_bs, vartype, loc){
  fish_result = fisher.test(table(data.frame(Disease = p_annotate_bs$VarType == vartype, Allosite = p_annotate_bs$allosite == loc)))
  fish_result
}
fish_bs(subset(allo_site, allosite %in% c("Binding Site", "Allo")), "Disease", "Allo")
fish_bs(subset(allo_site, allosite %in% c("Binding Site", "Allo")), "Polymorphism", "Allo")
fish_bs(subset(allo_site, allosite %in% c("Binding Site", "Allo")), "Unclassified", "Allo")

> length(unique(subset(allo_site, allosite %in% c("Binding Site", "Allo"))$UniProtID))
[1] 27

102 PDBs

> sum(!is.na(subset(allo_site,select = VarType)))
[1] 964


> with(allo_site, table(factor(proteinname)))                                                                                                                                                                    
"
Androgen receptor                                                                                                                                   
304                                                                                                                                   
Antithrombin-III                                                                                                                                   
436                                                                                                                                   
Coagulation factor IX                                                                                                                                   
307                                                                                                                                   
Cyclin-dependent kinase 2                                                                                                                                   
298                                                                                                                                   
Cytosolic purine 5'-nucleotidase                                                                                                                                   
470                                                                                                                                   
Dihydrofolate reductase                                                                                                                                   
186                                                                                                                                   
Dual specificity mitogen-activated protein kinase kinase 1                                                                                                                                   
316                                                                                                                                   
Farnesyl pyrophosphate synthase                                                                                                                                   
346                                                                                                                                   
Fructose-1,6-bisphosphatase 1                                                                                                                                   
319                                                                                                                                   
Glucokinase                                                                                                                                   
461                                                                                                                                   
Glutaminase kidney isoform, mitochondrial                                                                                                                                   
412                                                                                                                                   
Glycogen phosphorylase, liver form                                                                                                                                   
834                                                                                                                                   
Glycogen phosphorylase, muscle form                                                                                                                                   
821                                                                                                                                   
Hemoglobin subunit alpha                                                                                                                                   
191                                                                                                                                   
Hemoglobin subunit beta                                                                                                                                   
273                                                                                                                                   
Hexokinase-1                                                                                                                                   
903                                                                                                                                   
Integrin alpha-L                                                                                                                                   
184                                                                                                                                   
Kinesin-like protein KIF11                                                                                                                                   
351                                                                                                                                   
Leukotriene A-4 hydrolase                                                                                                                                   
610                                                                                                                                   
Mitogen-activated protein kinase 14                                                                                                                                   
359                                                                                                                                   
NAD-dependent malic enzyme, mitochondrial                                                                                                                                   
554                                                                                                                                   
Ribose-phosphate pyrophosphokinase 1                                                                                                                                   
309                                                                                                                                   
Serine/threonine-protein kinase Chk1                                                                                                                                   
279                                                                                                                                   
Transthyretin                                                                                                                                   
156                                                                                                                                   
Vitamin D3 receptor                                                                                                                                   
255                                                                                                                                   
[Pyruvate dehydrogenase (acetyl-transferring)] kinase isozyme 2, mitochondrial                                                                                                                                   
374                                                                                                                                   
cAMP and cAMP-inhibited cGMP 3',5'-cyclic phosphodiesterase 10A                                                                                                                                   
505 
"

[1] "HIS to ALA"
[1] "Disease"

## Bootstrap for confidence interval
library(boot)
aa_change_preference <- function(p_annotate_bs, interested_par){
  aa_change <- factor(unique(p_annotate_bs[,interested_par]))
  vartypes   = c("Disease", "Polymorphism", "Unclassified")
  ratio_amino <- data.frame(merge(aa_change, vartypes, all = TRUE))
  colnames(ratio_amino) <- c(interested_par, "vartypes")
  print(ratio_amino)
  for (vartype in vartypes) {
    for (each in aa_change){
      
      cannot_proc <- F
      
      oddsratio_fun <- function(d, i) {
        or <- fisher.test(table(data.frame(residueName = (d[i,interested_par] == each), Disease = d[i, "VarType"] == vartype)))$estimate[[1]]
        or
      }
      
      oddsboot <- tryCatch({
        boot(p_annotate_bs, oddsratio_fun, R = 1000, parallel = "multicore")
        
      }, error = function(err) { cannot_proc <<- T })
      
      oddsci   <- tryCatch({
        boot.ci(oddsboot, type = "norm")
      }, error = function(err) { cannot_proc <<- T })
      
      if (cannot_proc) {
        next
      }
      
      print(oddsci)
      
      ratio_amino[ratio_amino[,interested_par] == each & ratio_amino$vartypes == vartype,"estimate"] = oddsci$t0
      ratio_amino[ratio_amino[,interested_par] == each & ratio_amino$vartypes == vartype,"ci_low"] = oddsci$normal[2]
      ratio_amino[ratio_amino[,interested_par] == each & ratio_amino$vartypes == vartype,"ci_up"] = oddsci$normal[3]
    }
  }
  
  print(ratio_amino)
  # plot the graph
  ratio_amino_sort <- ratio_amino[with(ratio_amino, order(-estimate)),]
  ratio_amino[, interested_par] <- factor(ratio_amino[, interested_par], levels = as.character(subset(ratio_amino_sort, vartypes == "Disease")[, interested_par]))
  
  # for aa_change
  ratio_amino <- subset(ratio_amino, (ci_low - 1) * (ci_up - 1) > 0 & estimate > 0 & is.finite(estimate) & vartypes %in% c("Disease", "Polymorphism"))
  common_change <- intersect(subset(ratio_amino, vartypes == "Disease")$aa_change, subset(ratio_amino, vartypes == "Polymorphism")$aa_change)
  ratio_amino <- subset(ratio_amino, aa_change %in% common_change)
  
  ggplot(subset(ratio_amino, vartypes %in% c("Disease", "Polymorphism")), aes_string(x = interested_par, y = "estimate", ymin = "ci_low", ymax = "ci_up")) + geom_pointrange(aes(col = vartypes), position=position_dodge(width=0.30))  + ylab("Odds ratio & 95% CI") + geom_hline(aes(yintercept = 1)) + xlab("") + scale_y_log10()  
  ggsave(filename = "tmp.pdf", height=3, width=12)
  
  ratio_amino
}
ratio_amino <- aa_change_preference(protein_annotate_onlysnp, "uniprot_resnam_3d")

ratio_amino <- aa_change_preference(subset(protein_annotate_onlysnp, location %in% c("Binding Site", "Surface")), "uniprot_resnam_3d")

ratio_amino <- aa_change_preference(protein_annotate_onlysnp, "aa_change")

