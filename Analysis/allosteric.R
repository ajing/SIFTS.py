# allosteric
allo_annotate_org <- read.table("./Data/SNPOnStruct_biolip.txt", sep = "\t", header = F, quote = "", na.string = "\\N")

<<<<<<< HEAD
colnames(allo_annotate)<-c("pdbid", "biounit", "ModelID", "chainid", "resnam", "resnum", "structcode", "ssa", "rsa", "UniProtID", "uniprot_resnam", "uniprot_resnum", "bs_biounit", "bs_p_chainid", "bs_p_resnum", "ligandName", "allo_pubmedId", "allo_pubmedTitle", "allo_organism", "allo_bsid", "genename", "SwissProt_AC", "FTID", "AABefore", "prot_resnum", "AAAfter", "VarType", "dbSNPID", "DiseaseName", "uniprot_resnam_3d", "AABeforeProp", "AAAfterProp", "proteinname", "reviewed", "gene_name_acc", "interproname_acc")
=======
colnames(allo_annotate_org)<-c("pdbid", "biounit", "ModelID", "chainid", "resnam", "resnum", "structcode", "ssa", "rsa", "UniProtID", "uniprot_resnam", "uniprot_resnum", "bs_chainid", "bs_resnum", "ligandName", "PubMedID", "PubMedTitle", "Organism", "BSID", "genename", "SwissProt_AC", "FTID", "AABefore", "prot_resnum", "AAAfter", "VarType", "dbSNPID", "DiseaseName", "uniprot_resnam_3d", "AABeforeProp", "AAAfterProp", "proteinname", "reviewed", "gene_name_acc", "interproname_acc")

library(sqldf)
allo_annotate = sqldf("select *, max(ssa) as ssa_m from allo_annotate_org group by UniProtID, uniprot_resnum, dbSNPID")
>>>>>>> d4408a32e0c41ef84d2725db219051bc0ef4ee17

allo_annotate$location = "Core"

#MAX_ACC["GLY"] = 84.0
<<<<<<< HEAD
allo_annotate$location[allo_annotate$ssa > 5] = "Surface"

allo_annotate[ which(allo_annotate$distance <= 4.0), 'location'] = "Binding Site"

allo_annotate[ which(!is.na(allo_annotate$allo_pubmedId)), 'location'] = "Allo"
# For allosteric site
allo_annotate[]
=======
allo_annotate$location[allo_annotate$rsa > 0.05] = "Surface"
>>>>>>> d4408a32e0c41ef84d2725db219051bc0ef4ee17

allo_annotate[!is.na(allo_annotate$ligandName), "location"] = "Binding Site"

# only keep UniProt with at least one nsSNP
require(plyr)
#snp_num = ddply(allo_annotate, .(UniProtID), summarise, N = length(unique(uniprot_resnum[!is.na(VarType)])))
allo_annotate_withsnp = subset(allo_annotate, UniProtID %in% factor(unique(subset(allo_annotate, !is.na(VarType))$UniProtID)))
length(unique(allo_annotate_withsnp$UniProtID))
"42"

table(allo_annotate_withsnp$VarType, allo_annotate_withsnp$location)

"Binding Site Core Surface
Disease              1739  125     444
Polymorphism          288    8     129
Unclassified          682   10     367"

table(allo_annotate_withsnp$location)


length(unique(allo_annotate[allo_annotate$VarType == "Unclassified",]$UniProtID)) #17

length(unique(allo_annotate[allo_annotate$VarType == "Disease",]$UniProtID))  #24

length(unique(allo_annotate[allo_annotate$VarType == "Polymorphism",]$UniProtID))  #34

# allosteric sites or other sites
<<<<<<< HEAD
=======
allo_annotate_withsnp$allosite = allo_annotate_withsnp$location
allo_annotate_withsnp$allosite[!is.na(allo_annotate_withsnp$PubMedID)] = "Allo"

table(allo_annotate_withsnp$VarType, allo_annotate_withsnp$allosite)

"Allo Binding Site Core Surface
Disease       108         1633  125     442
Polymorphism    6          288    8     123
Unclassified   64          638   10     347"

subset( allo_annotate_withsnp, select = c(resnum, location, VarType, ssa, rsa, ligandName), UniProtID == "O14757")

allo_mix = rbind(subset( allo_site, select = c(UniProtID, uniprot_resnum, location, VarType, allosite)), subset( allo_annotate_withsnp, select = c(UniProtID, uniprot_resnum, location, VarType, allosite), !(UniProtID %in% allo_site$UniProtID)))


sort(allo_annotate_withsnp$UniProtID)
tmp = merge( subset( allo_site, select = c(resnum, location, VarType), UniProtID == "O14757"), subset( allo_annotate_withsnp, select = c(resnum, location, VarType), UniProtID == "O14757"), by = "resnum")

tmp = merge( subset( allo_site, select = c(UniProtID, uniprot_resnum, location, VarType, allosite)), subset( allo_annotate_withsnp, select = c(UniProtID, uniprot_resnum, location, VarType, allosite)), by = c("uniprot_resnum", "UniProtID"))

tmp = merge( subset( allo_site, select = c(UniProtID, uniprot_resnum, location, VarType, allosite)), subset( allo_mix, select = c(UniProtID, uniprot_resnum, location, VarType, allosite)), by = c("uniprot_resnum", "UniProtID"))
table(tmp$location.x, tmp$location.y)
table(allosite = tmp$allosite.x, newdata = tmp$allosite.y)


subset(tmp, location.y == "Binding Site", select = c(location.x, location.y))
summary(subset(tmp, location.y == "Binding Site", select = c(location.x, location.y)))

>>>>>>> d4408a32e0c41ef84d2725db219051bc0ef4ee17
fish_bs <- function(p_annotate_bs, vartype, loc){
  fish_result = fisher.test(table(data.frame(Disease = p_annotate_bs$VarType == vartype, Allosite = p_annotate_bs$allosite == loc)))
  fish_result
}
fish_bs(subset(allo_annotate_withsnp, allosite %in% c("Binding Site", "Allo")), "Disease", "Allo")

fish_bs(subset(allo_mix, allosite %in% c("Binding Site", "Allo")), "Disease", "Allo")


"
        Fisher's Exact Test for Count Data

data:  
p-value = 0.4207
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.6343957 1.2128784
sample estimates:
odds ratio 
 0.8749319"

fish_bs(subset(allo_site, allosite %in% c("Binding Site", "Allo")), "Polymorphism", "Allo")

fish_bs(subset(allo_mix, allosite %in% c("Binding Site", "Allo")), "Polymorphism", "Allo")

"        Fisher's Exact Test for Count Data

data:  
p-value = 0.0003824
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.09868007 0.61940961
sample estimates:
odds ratio 
  0.275155 "

fish_bs(subset(allo_site, allosite %in% c("Binding Site", "Allo")), "Unclassified", "Allo")

<<<<<<< HEAD

# comparing allosteric sites and other ligand binding sites for each protein
each_uniprot = by(subset(allo_site, allosite %in% c("Binding Site", "Allo")), subset(allo_site, allosite %in% c("Binding Site", "Allo"))[, "UniProtID"], function(x) {fish_bs(x, "Disease", "Allo")})


snp_num = ddply(allo_site, .(UniProtID), summarise, fish_odds = fisher.test(table(data.frame(Disease = VarType, Allosite = allosite)))$estimate, fish_pvalue = fisher.test(table(data.frame(Disease = VarType, Allosite = allosite)))$p.value)
=======
"        Fisher's Exact Test for Count Data

data:  
p-value = 0.001797
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 1.208083 2.348072
sample estimates:
odds ratio 
  1.690005 "
>>>>>>> d4408a32e0c41ef84d2725db219051bc0ef4ee17
