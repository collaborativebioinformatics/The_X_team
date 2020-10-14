library(stringr)

args = commandArgs(trailingOnly=TRUE)

t <- read.table(args[1])

head(t)

hom <- which(str_detect(t$V1, "Hom"))
het <- which(str_detect(t$V1, "Het"))
low <- which(str_detect(t$V1, "Lowfreq"))
maxLen <- 1e4
shortSV <- which(t$V9 < maxLen & t$V9 >= 50)
matConf <- which(t$V7 < 10 | t$V7 < 0.05*t$V9)
patConf <- which(t$V15 < 10 | t$V15 < 0.05*t$V17)

matConfStrict <- which(t$V7 < 5 | t$V7 < 0.05*t$V9)
patConfStrict <- which(t$V15 < 5 | t$V15 < 0.05*t$V17)




both <- union(matConf,patConf)
notConf <- setdiff(seq(1,length(t$V15)), both)
br <- seq(0,maxLen,100)
h1 <- hist(t$V9[intersect(shortSV,intersect(hom,both))],breaks=br, plot=F)
h2 <- hist(t$V9[intersect(shortSV,intersect(het,both))],breaks=br, plot=F)
h3 <- hist(t$V9[intersect(shortSV,intersect(low,both))],breaks=br, plot=F)

hnc1 <- hist(t$V9[intersect(shortSV,intersect(hom,notConf))],breaks=br, plot=F)
hnc2 <- hist(t$V9[intersect(shortSV,intersect(het,notConf))],breaks=br, plot=F)
hnc3 <- hist(t$V9[intersect(shortSV,intersect(low,notConf))],breaks=br, plot=F)



yr <- range(c(h1$counts,h2$counts,h3$counts,hnc1$counts,hnc2$counts,hnc3$counts))
yr <- c(1,max(yr))
library(RColorBrewer)
p <- brewer.pal(3,"Set1")
pdf("CombinedHist.pdf")
plot(h1$mids, h1$counts+1, ylim=yr, type='l', lty=1, col=p[1], log="y")
points(h2$mids, h2$counts, ylim=yr, type='l', lty=1, col=p[2])
points(h3$mids, h3$counts, ylim=yr, type='l', lty=1, col=p[3])

points(hnc1$mids, hnc1$counts, ylim=yr, type='l', lty=2, col=p[3])
points(hnc2$mids, hnc2$counts, ylim=yr, type='l', lty=2, col=p[2])
points(hnc3$mids, hnc3$counts, ylim=yr, type='l', lty=2, col=p[3])
legend("topright", legend=c("Homozygous validated", "Heterozygous validated", "Low-cov validated", "Homozygous not validated", "Heterozygous not validated", "Low-cov not validated"), lty=c(1,1,1,2,2,2), col=p)
dev.off()
