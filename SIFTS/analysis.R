filename = "../pair.txt"
#filename = "pair.txt"
xy = read.table(filename, sep = "\t")
x  = xy[,1]
y  = xy[,2]

install.packages("hexbin")

library(hexbin)
hbo <- hexbin(x,y, xlab = "spatial distance", ylab = "sequence distance")
plot(hbo)
jpeg("hbo.jpg")
plot(hbo)
dev.off()

hist(x)
histpercent("pairdist.jpg", x, "Pair Spatial Distance", "Distance in Angstrom", "Probability")

content <- data.frame
content <- read.table("order_join.txt", sep="\t", quote = "")
pie(content[,2], labels = content[,1])

distfile <- "../onlydist.txt_noinf"
distance <- read.table(distfile, header=F)

distance_dat <- data.matrix(distance)
distance_num <- as.numeric(distance[,"V1"])


histpercent <- function(filename, datavector, gtitle, xlabel, ylabel){
jpeg(filename)
disthist = hist(datavector)
disthist$counts = disthist$counts / sum(disthist$count)
plot(disthist, main = gtitle, ylab = ylabel, xlab = xlabel)
dev.off()
}

histpercent("bsdist_hist.jpg", distance[distance[,1]<250,1], "Distance to Binding Site", "Distance in Angstrom", "Probability")
histpercent("bsdist_hist.jpg", distance[,1], "Distance to Binding Site",  "Distance in Angstrom", "Probability")

# For distance between residues
filename = "../distaa.txt_filtered"
distdata = read.table(filename)

histpercent("bsdist_hist_allchains_snp.jpg", as.numeric(distdata[,"V7"]), "Distance to Binding Site", "Distance in Angstrom", "Probability")
histpercent("bsdist_hist_allchains_less10_snp.jpg", as.numeric(distdata[distdata[,"V7"]< 10,"V7"]), "Distance to Binding Site", "Distance in Angstrom", "Probability")
histpercent("bsdist_hist_allchains_less20_snp.jpg", as.numeric(distdata[distdata[,"V7"]< 20,"V7"]), "Distance to Binding Site", "Distance in Angstrom", "Probability")

filename = "../distaa.txt"
distdata = read.table(filename)
histpercent("bsdist_hist_allchains.jpg", as.numeric(distdata[,"V7"]), "Distance to Binding Site", "Distance in Angstrom", "Probability")
histpercent("bsdist_hist_allchains_less10.jpg", as.numeric(distdata[distdata[,"V7"]< 10,"V7"]), "Distance to Binding Site", "Distance in Angstrom", "Probability")
histpercent("bsdist_hist_allchains_less20.jpg", as.numeric(distdata[distdata[,"V7"]< 20,"V7"]), "Distance to Binding Site", "Distance in Angstrom", "Probability")


disease_snp <- read.table("../disease_snp_dist.txt", header=F)
disease_snp <- data.matrix(disease_snp)

datavector <- disease_snp
plot(disthist,  xlim=c(0,200))

histpercent("disease_dist.jpg", disease_snp, "Distance from Disease SNP to BindingSite", "Probability", "Distance in Angstrom")
histpercent("disease_dist_less10.jpg", disease_snp[disease_snp < 10], "Distance from Disease SNP to BindingSite", "Probability", "Distance in Angstrom")

library(ggplot2)
disease_snp <- read.table("../disease_snp_dist.txt", header=F)
disease_snp <- as.numeric(unlist(disease_snp))

snp_dist = data.frame(distance = as.numeric(distdata[,"V7"]))
disease_dist = data.frame(distance = disease_snp)
snp_dist$type = "all_snp"
disease_dist$type = "disease_snp"
ggplot(rbind(snp_dist,disease_dist), aes(distance , fill = type)) + geom_density(alpha = 0.2)
