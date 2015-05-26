# dual (or triple) mutational distribution from COSMIC

# The dual mutataion distribution on UniProt disease nsSNP and polymorphism nsSNP

## load data from humsavar.txt
humsavar <- read.table("/users/ajing/SNPDist/Data/humsavar_tru.txt.out_test", sep = "\t", header = F, quote = "", na.string = "-")
humsavar <- read.table("humsavar_tru.txt.out_test", sep = "\t", header = F, quote = "", na.string = "-")


colnames(humsavar) <- c("genename", "uniprot_id", "var_id", "aa_before", "resnum", "aa_after", "var_type", "snp_id", "disease")

dist_each_uniprot <- merge(humsavar, humsavar, by = "uniprot_id")
dist_each_uniprot$dist <- dist_each_uniprot$resnum.x - dist_each_uniprot$resnum.y
dist_each_uniprot <- subset(dist_each_uniprot, dist > 0)
dist_each_uniprot$var_type_c <- apply(dist_each_uniprot, 1, function(x){paste(sort(c(x[["var_type.x"]],x[["var_type.y"]])), collapse = " to ")})

# The mutation distribution between disease nsSNP and polymorphism nsSNP
require(ggplot2)
ggplot(subset(dist_each_uniprot, var_type_c %in% c("Disease to Disease", "Polymorphism to Polymorphism", "Unclassified to Unclassified")), aes(x = dist, color = var_type_c)) + scale_x_log10() + geom_density()
ggsave("snp_dist.pdf")

# The mutation distribution between cancer disease nsSNP and polymorphism nsSNP
ggplot(rbind(subset(dist_each_uniprot, var_type_c %in% c("Disease to Disease", "Polymorphism to Polymorphism", "Unclassified to Unclassified") & grepl("cancer", disease.x) & grepl("cancer", disease.y)), subset(dist_each_uniprot, var_type_c %in% c("Polymorphism to Polymorphism"))), aes(x = dist, color = var_type_c) ) + scale_x_log10() + geom_density()
ggsave("snp_dist_cancer.pdf")
