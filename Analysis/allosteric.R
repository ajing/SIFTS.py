# allosteric
allo_annotate<- read.table("./Data/SNPOnStruct_biolip.txt", sep = "\t", header = F, quote = "", na.string = "\\N")

colnames(allo_annotate)<-c("pdbid", "biounit", "ModelID", "chainid", "resnam", "resnum", "structcode", "ssa", "rsa", "UniProtID", "uniprot_resnam", "uniprot_resnum", "bs_chainid", "bs_resnum", "ligandName", "PubMedID", "PubMedTitle", "Organism", "BSID", "genename", "SwissProt_AC", "FTID", "AABefore", "prot_resnum", "AAAfter", "VarType", "dbSNPID", "DiseaseName", "uniprot_resnam_3d", "AABeforeProp", "AAAfterProp", "proteinname", "reviewed", "gene_name_acc", "interproname_acc")

allo_annotate$location = "Core"

#MAX_ACC["GLY"] = 84.0
allo_annotate$location[allo_annotate$ssa > 5] = "Surface"

allo_annotate[!is.na(allo_annotate$ligandName), "location"] = "Binding Site"

table(allo_annotate$VarType, allo_annotate$location)

Binding Site Core Surface
Disease              1739  125     444
Polymorphism          288    8     129
Unclassified          682   10     367


# only keep UniProt with at least one nsSNP
require(plyr)
#snp_num = ddply(allo_annotate, .(UniProtID), summarise, N = length(unique(uniprot_resnum[!is.na(VarType)])))
allo_annotate_withsnp = subset(allo_annotate, UniProtID %in% factor(unique(subset(allo_annotate, !is.na(VarType))$UniProtID)))
length(unique(allo_annotate_withsnp$UniProtID))

table(allo_annotate_withsnp$VarType, allo_annotate_withsnp$location)

"Binding Site Core Surface
Disease              1739  125     444
Polymorphism          288    8     129
Unclassified          682   10     367"


length(unique(allo_annotate[allo_annotate$VarType == "Unclassified",]$UniProtID)) #17

length(unique(allo_annotate[allo_annotate$VarType == "Disease",]$UniProtID))  #24

length(unique(allo_annotate[allo_annotate$VarType == "Polymorphism",]$UniProtID))  #34

# allosteric sites or other sites
allo_res <- read.table("./Data/pro_res_count.txt", sep = "\t", header = F, quote = "", na.string = "\\N")
colnames(allo_res) <- c("UniProtID", "uniprot_resnum")
allo_res$allosite = "Allo"
allo_site <- merge(allo_annotate_withsnp, allo_res, all.x = T)
allo_site[is.na(allo_site$allosite), "allosite"] <- allo_site[is.na(allo_site$allosite), "location"]
allo_site$allosite <- factor(allo_site$allosite)

table(allo_site$VarType, allo_site$allosite)

"Allo Binding Site Core Surface
Disease       108         1633  125     442
Polymorphism    6          288    8     123
Unclassified   64          638   10     347"



fish_bs <- function(p_annotate_bs, vartype, loc){
  fish_result = fisher.test(table(data.frame(Disease = p_annotate_bs$VarType == vartype, Allosite = p_annotate_bs$allosite == loc)))
  fish_result
}
fish_bs(subset(allo_site, allosite %in% c("Binding Site", "Allo")), "Disease", "Allo")

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

"        Fisher's Exact Test for Count Data

data:  
p-value = 0.001797
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 1.208083 2.348072
sample estimates:
odds ratio 
  1.690005 "