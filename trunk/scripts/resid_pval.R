genes_pvalues <- read.delim("../data/genes_pvalues.txt", header=F)
# BF_pval <- ''
 BF_pval <- merge(exon_BF_bins, genes_pvalues, by.x = "Bin", by.y = "V1")
 View(BF_pval)

# for (i in 1:length(exon_BF_bins)){
#   if (exon_BF_bins[i,5] >= .2){
#     exon_BF_bins[i,5] <- -1
#   }
# }

cor.p <- cor.test(BF_pval$resid,BF_pval$V2)
plot(BF_pval$resid,
     axes=FALSE, ann=FALSE,
     ylab = "Bin residuals", 
     xlab = "Bins",
     xlim=c(0,262), 
     ylim=c(0,0.2),
     pch = 20, 
     cex = .3, 
     col = "blue")
par(new=TRUE)
plot(BF_pval$V2, col = "green",
     axes=FALSE, ann=FALSE,
     ylab = "Bin residuals", 
     xlab = "Bins",
     xlim=c(0,262), 
     ylim=c(0,0.2))
axis(1, las=2, at=0:(length(BF_pval$Bin)-1), lab=BF_pval$Bin, cex.axis=.15)
axis(2)
box()
# Graph resid with thicker blue dashed line
lines(BF_pval$resid, type="l", lty=2, lwd=2, 
  col="blue")

# Graph pvalue with thicker green dotted line
lines(BF_pval$V2, type="l", lty=3, lwd=2, 
  col="green")
title(main = "Residuals for each bin")