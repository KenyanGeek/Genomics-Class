setwd("/Users/chege/Desktop/Spring 2020/SCHOOL/Genomics Class/thetas")

list.files()


SFS <- scan("XCS_outfold.sfs")

sumSFS <- sum(SFS)

pctPoly = 100*(1-(SFS[1]/sumSFS))

plotSFS <- SFS[-c(1,length(SFS))]             
                
barplot(plotSFS)

div <- read.table("XCS_folded_allsites.thetas.idx.pestPG")

colnames(div) = c("window","chrname","wincenter","tW","tP","tF","tH","tL","tajD","fulif","fuliD","fayH","zengsE","numSites")

div$tWpersite = div$tW/div$numSites
div$tPpersite = div$tP/div$numSites


pdf("XCS_diversity_stats.pdf")
par(mfrow=c(2,2))
hist(div$tWpersite, col = "gray", xlab = "Theta-W", main = "")
hist(div$tPpersite, col = "gray", xlab = "Theta-P", main = "")
hist(div$tajD, col = "gray", xlab = "Tajima's D", main = "")
barplot(plotSFS)
dev.off()

summary(div)

Tajd <- div$tPpersite - div$tWpersite

summary(Tajd)
