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
sc_volumn$volumn.c = cut(sc_volumn$volumn, breaks=3, label = c("small", "median", "large"))
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