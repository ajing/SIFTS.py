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
filename = "../distaa.txt"
distdata = read.table(filename)
all_dist = data.frame(distance = as.numeric(distdata[,"V7"]))

disease_snp <- read.table("../disease_snp_dist.txt", header=F)
disease_snp <- as.numeric(unlist(disease_snp))

polymore_snp <- read.table("/tmp/polymorephism_snp_dist.txt", header=F)
polymore_snp <- as.numeric(unlist(polymore_snp))

filename = "../distaa.txt_filtered"
distdata = read.table(filename)
snp_dist = data.frame(distance = as.numeric(distdata[,"V7"]))
disease_dist = data.frame(distance = disease_snp)
snp_dist$type = "all_snp"
disease_dist$type = "disease_snp"


all_dist$type = "all_pair"
jpeg("all_hist.jpg")
ggplot(rbind(all_dist, snp_dist,disease_dist), aes(distance , fill = type)) + geom_density(alpha = 0.2)
dev.off()


polymore_dist = data.frame(distance = polymore_snp)
polymore_dist$type = "polymorephism_snp"
jpeg("polymore_hist.jpg")
ggplot(rbind(snp_dist,polymore_dist), aes(distance , fill = type)) + geom_density(alpha = 0.2)
dev.off()

unclass_snp <- read.table("/tmp/unclassified_snp_dist.txt", header=F)
unclass_snp <- as.numeric(unlist(unclass_snp))
unclass_dist = data.frame(distance = unclass_snp)
unclass_dist$type = "unclassified_snp"
jpeg("unclass_hist.jpg")
ggplot(rbind(snp_dist,unclass_dist), aes(distance , fill = type)) + geom_density(alpha = 0.2)
dev.off()


# The first column is family size (here is the number of pairs.. so very big number) second is max sequence length third column is standard deviation of distance
dist_var <- read.table("../Data/dist_variance.txt", sep="\t", quote= "", header=F)
plot(dist_var[,1], dist_var[,3], xlab = "family size", ylab = "standard deviation of distance", col = dist_var[,2])

install.packages("plot3D")
library(plot3D)
scatter3D(dist_var[,1], dist_var[,2], dist_var[,3])
pdf("dist_sd_family_relation.pdf")
png("dist_sd_family_relation.png")
scatter2D(dist_var[,1], dist_var[,3], colvar=dist_var[,2], log = "x", xlab = "The number of distances for a family (to binding sites)", ylab = "standard deviation of distances", clab = "Sequence Length")
index = dist_var[,3] > 60 & dist_var[,1] / dist_var[,2] ** 2 > 0.1
text(dist_var[index,1], dist_var[index,3], dist_var[index,4], cex=0.7, pos=1, col="black")
dev.off()
