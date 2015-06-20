# To convey the idea that mutagenesis and existing methods are different

# the data is from PredictSNP
predict_snp <- read.table("../Data/pathogenicity.csv", header = T, sep = "\t", quote = "")

predict_snp$Target.residue <- ordered(predict_snp$Target.residue, levels = rev(levels(predict_snp$Target.residue)))

library(ggplot2)
ggplot(predict_snp, aes(Position, Target.residue, fill = PredictSNP.confidence)) + geom_raster() + scale_x_continuous(breaks = round(seq(min(predict_snp$Position), max(predict_snp$Position), by = 1),1)) + scale_fill_gradient(low = "white", high = "steelblue") + theme_bw()

ggsave("mutagenesis.jpg", width = 15, height = 4)
