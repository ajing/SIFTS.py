library(ggplot2)
library(plyr)

surface <- read.table("./Result/avedist2surface.txt_small", sep = "\t", header = F)

colnames(surface)<-c("pdbid", "chainid", "resnum", "resnam", "surf_dist", "surf_dist_alpha")

surface <- melt(surface, id = c("pdbid","chainid","resnam","surf_dist","sur_dist_alpha"))
library("reshape2")
surface <- melt(surface, id = c("pdbid","chainid","resnum", "resnam"))
qplot(value, data = surface, color = variable, geom="density", log="x")
newsurf<-subset(surface, value != "None")
subset(newsurf, is.na(as.numeric(value)))
surface <- newsurf


qplot(value, data = surface, geom="density", log="x") + facet_wrap(resnam~.)
aanam <- c('ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLU', 'GLN', 'GLX', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL')
qplot(value, data = subset(surface,resnam %in% aanam[1:11]), color = resnam, geom="density") + coord_cartesian(xlim = c(1,10))
qplot(value, data = subset(surface,resnam %in% aanam[12:22]), color = resnam, geom="density")  + coord_trans(x="log2")+ coord_cartesian(xlim = c(1,10))


library("scales")
qplot(value, data = subset(surface,resnam %in% aanam[12:22]), color = resnam, geom="density")  + scale_x_continuous(trans=log2_trans())+ coord_cartesian(xlim = c(1,10))
snp_dist <-read.csv("./Result/snp_dist.txt",sep="\t", header=F)
colnames(snp_dist)<-c("pdbid","snpid1","chainid1","resnam1","resnum1","snpid2","chainid2","resnam2","resnum2", "spa_dist","seq_dist")

snp_dist_nofalse <- subset(snp_dist, seq_dist != 'False')
snp_dist_nofalse[,"seq_dist"] <- as.numeric(levels(snp_dist_nofalse$seq_dist)[snp_dist_nofalse$seq_dist])
ggplot(snp_dist_nofalse, aes(x=spa_dist, y=seq_dist)) + geom_bin2d()

humsavar <- read.table("./Data/humsavar_tru.txt.out", sep="\t",header=F, quote="")
humsavar <- subset(humsavar, V6 !="-")
humsavar <- humsavar[,c("V5","V6")]
colnames(humsavar)<-c("type", "snpid")

snp_dist_nof_var<-merge(x = snp_dist_nofalse, y = humsavar, by.x = "snpid1", by.y="snpid", all.x=TRUE)
colnames(snp_dist_nof_var)[12]<- "type1"
snp_dist_nof_var<-merge(x = snp_dist_nof_var, y = humsavar, by.x = "snpid2", by.y="snpid", all.x=TRUE)
colnames(snp_dist_nof_var)[13]<- "type2"
ggplot(subset(snp_dist_nof_var,type1 != "Unclassified" & type2 != "unclassified"), aes(x=spa_dist, y=seq_dist)) + geom_bin2d() + facet_grid(type1 ~.)
ggplot(subset(snp_dist_nof_var,type1 != "Unclassified" & type2 != "Unclassified"), aes(x=spa_dist, y=seq_dist)) + stat_density2d(geom="tile",aes(fill=..density..), contour=FALSE) + facet_grid(type1 + type2 ~.)

snp_dist_nof_var$type = paste(snp_dist_nof_var$type1,snp_dist_nof_var$type2)
snp_dist_nof_var[snp_dist_nof_var$type1 != snp_dist_nof_var$type2 & snp_dist_nof_var$type1 != "Unclassified" & snp_dist_nof_var$type2 != "Unclassified","type"] <- "Polymorphism Disease"
snp_dist_nof_var[snp_dist_nof_var$type1 != snp_dist_nof_var$type2 & snp_dist_nof_var$type1 != "Polymorphism" & snp_dist_nof_var$type2 != "Polymorphism","type"] <- "Unclassified Disease"
snp_dist_nof_var[snp_dist_nof_var$type1 != snp_dist_nof_var$type2 & snp_dist_nof_var$type1 != "Disease" & snp_dist_nof_var$type2 != "Disease","type"] <- "Polymorphism Unclassified"
snp_dist_nof_var$type = as.factor(snp_dist_nof_var$type)

ggplot(subset(snp_dist_nof_var,type1 != "Unclassified" & type2 != "Unclassified"), aes(x=spa_dist, y=seq_dist)) + stat_density2d(geom="tile",aes(fill=..density..), contour=FALSE) +  scale_fill_gradient(low = "white", high = "blue") + facet_grid(type ~.) + coord_cartesian(xlim = c(0, 50))
ggplot(subset(snp_dist_nof_var,type1 != "Unclassified" & type2 != "Unclassified"), aes(x=spa_dist, y=seq_dist))  + facet_grid(type ~.) + stat_density2d(geom="tile",aes(fill=..density..), contour=FALSE) +  scale_fill_gradient(low = "white", high = "blue")+ scale_x_continuous(limits = c(0, 50), expand=c(0,0)) + scale_y_continuous(limits = c(0,200), expand=c(0,0))

ggplot(subset(snp_dist_nof_var,type1 != "Unclassified" & type2 != "Unclassified" & type %in% c("Polymorphism Polymorphism")), aes(x=spa_dist, y=seq_dist))  + facet_grid(type ~.) + stat_density2d(geom="tile",aes(fill=..density..), contour=FALSE) +  scale_fill_gradient(low = "white", high = "blue")+ scale_x_continuous(limits = c(0, 50), expand=c(0,0)) + scale_y_continuous(limits=c(0,200), expand=c(0,0))

dist_bs <- read.table("./Result/distaa_biolip.txt_filtered", sep = "\t", header =F)
colnames(dist_bs)<-c("pdbid","bs_chainid","bsid","resnam","res_chainid","resnum", "distance", "snpid")
dist_bs$pdbid<-toupper(dist_bs$pdbid)

disease_spot<-subset(snp_dist_nof_var,type1 != "Unclassified" & type2 != "Unclassified" & type %in% c("Disease Disease"))
disease_spot[!(disease_spot$spa_dist < 17 & disease_spot$spa_dist > 7 & disease_spot$seq_dist < 20), "spot"] <- "outside"
disease_spot[disease_spot$spa_dist < 13 & disease_spot$spa_dist > 7 & disease_spot$seq_dist < 20, "spot"] <- "spot1"
disease_spot[disease_spot$spa_dist < 17 & disease_spot$spa_dist >= 13 & disease_spot$seq_dist < 20, "spot"] <- "spot2"
disease_spot1<- disease_spot[,c("pdbid", "snpid1", "chainid1","resnam1","resnum1","spot")]
disease_spot2<- disease_spot[,c("pdbid", "snpid2", "chainid2","resnam2","resnum2","spot")]
colnames(disease_spot1)<- c("pdbid", "snpid", "chainid","resnam","resnum","spot")
colnames(disease_spot2)<- c("pdbid", "snpid", "chainid","resnam","resnum","spot")
disease_spot_melt<-rbind(disease_spot1,disease_spot2)

disease_spot_melt_dist<-merge(x = disease_spot_melt, y = dist_bs[,c("pdbid","snpid","res_chainid", "resnum","distance")], by.x = c("pdbid","snpid","chainid","resnum"), by.y = c("pdbid","snpid","res_chainid","resnum"), all.x=TRUE)


c_disease_spot <- ddply(subset(disease_spot_melt_dist, !is.na(distance)), "spot", summarise, mean=mean(distance))
qplot(distance, data=subset(disease_spot_melt_dist, !is.na(distance)), geom="density", facets=spot~., alpha=I(0.3)) + geom_vline(data=c_disease_spot, aes(xintercept=mean,  colour=spot),linetype="dashed", size=1)+ xlim(c(0,100))



pdf("eachPDB.pdf")
ggplot(subset(snp_dist_nof_var,type1 != "Unclassified" & type2 != "Unclassified"), aes(x=spa_dist, y=seq_dist))  + facet_grid(type ~ pdbid) + stat_density2d(geom="tile",aes(fill=..density..), contour=FALSE) +  scale_fill_gradient(low = "white", high = "blue")+ scale_x_continuous(limits = c(0, 50), expand=c(0,0)) + scale_y_continuous(limits=c(0,200), expand=c(0,0))
dev.off()


pdf("eachpdb.pdf")
for (pdb in sample(unique(subset(snp_dist_nof_var, type1 != "Unclassified" & type2 != "Unclassified")$pdbid),10)) {
  print(dim(subset(snp_dist_nof_var,type1 != "Unclassified" & type2 != "Unclassified" & pdbid == pdb)))
  pdbplot <- ggplot(subset(snp_dist_nof_var,type1 != "Unclassified" & type2 != "Unclassified" & pdbid == pdb), aes(x=spa_dist, y=seq_dist, color=type))  + geom_point() 
  print(pdbplot)
}
dev.off()

split_pdb <- function(d, tot_length) {
  x <- seq_along(d)
  split(d, ceiling(x/tot_length))
}

pdbplot <- ggplot(subset(snp_dist_nof_var,type1 != "Unclassified" & type2 != "Unclassified"), aes(x=spa_dist, y=seq_dist, color=type))  + facet_wrap(~pdbid) + geom_point()
ggsave(plot=pdbplot, filename="all_pdb.pdf", height=49, width=49)



pdb_col <- split_pdb(unique(snp_dist_nof_var$pdbid), 200)
for (index_col in names(pdb_col)) {
  print(dim(subset(snp_dist_nof_var,type1 != "Unclassified" & type2 != "Unclassified" & pdbid %in% lapply(pdb_col[index_col],as.character)[[index_col]])))
  pdbplot <- ggplot(subset(snp_dist_nof_var,type1 != "Unclassified" & type2 != "Unclassified" & pdbid %in% lapply(pdb_col[index_col],as.character)[[index_col]]), aes(x=spa_dist, y=seq_dist, color=type)) + facet_wrap(~pdbid) + geom_point() 
  ggsave(plot=pdbplot, filename=paste(index_col, ".jpeg", sep=""), height=20, width=25)
}


ggplot(subset(snp_dist_nof_var,type1 != "Unclassified" & type2 != "Unclassified" & pdbid == "1E2S"), aes(x=spa_dist, y=seq_dist))  + facet_grid(type ~.) + stat_density2d(geom="tile",aes(fill=..density..), contour=FALSE) +  scale_fill_gradient(low = "white", high = "blue")+ scale_x_continuous(limits = c(0, 50), expand=c(0,0)) + scale_y_continuous(limits=c(0,200), expand=c(0,0))


convert_prop <-read.csv("./Result/convert_property_snp.txt",sep="\t", header=F)
colnames(convert_prop) <- c("type1", "type2","count","type","snp_type")
sum_convert<-aggregate(count~type2 + snp_type, data = subset(convert_prop, type != "\\N"), FUN = sum)
ggplot(data=subset(sum_convert, snp_type %in% c("Disease", "Polymorphism")), aes(x=factor(1), y=count, fill=factor(type2))) + geom_bar(stat="identity") + facet_grid(snp_type~., scale="free_y") + coord_polar(theta="y")


snp_annotate <-read.csv("./Data/snpannotate.txt",sep="\t", header=F)
colnames(snp_annotate) <- c("pdbid", "snpid", "chainid","snp_position","ref_residue","alt_residue", "chr", "strand", "chr_position", "secondary", "solubility", "ligandcodes", "conservation", "domain_interface", "genename", "swissprot_ac", "ftid", "AAChange", "snptype", "dbsnpid", "diseasename")
qplot(conservation, data = subset(snp_annotate, snptype %in% c("Disease", "Polymorphism")) , color = snptype, geom = "density")


ggplot(data=subset(snp_annotate, snptype %in% c("Disease", "Polymorphism")), aes(x=factor(1), y=chr, fill=factor(chr))) + geom_bar() 
