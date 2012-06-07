exon_BF_bins <- read.csv("../data/combined_output_b37-BinFreq.csv")

resid <- abs(exon_BF_bins$Control.Freq. - exon_BF_bins$Case.Freq.)/sqrt(2)
exon_BF_bins <- cbind(exon_BF_bins, resid)
for (i in 1:length(exon_BF_bins$Bin)){
  if (exon_BF_bins[i,2]<0 || exon_BF_bins[i,3]<0){
    exon_BF_bins[i,4] <- -1
  }
  if (exon_BF_bins[i,1] == "ATP8B4" || exon_BF_bins[i,1] =="DUSP26" || exon_BF_bins[i,1] =="ARHGEF10L"){
    show(exon_BF_bins[i,])
  }
}
View(exon_BF_bins)

#find the limit of the highest ten residuals 
sorted_resid <- sort(exon_BF_bins$resid, TRUE)
resid_lim <- sorted_resid[30]

#create a column with ids for only those top 10 residuals
labels_top <- exon_BF_bins$Bin
top10table <- exon_BF_bins[0,]
for (i in 1:length(exon_BF_bins$Bin)){
  if (exon_BF_bins[i,4] < resid_lim){
    labels_top[i] <- ''
  }
  else{
    top10table <- rbind(top10table,exon_BF_bins[i,])
  }
}
show(top10table)
    
plot(exon_BF_bins$Control.Freq., exon_BF_bins$Case.Freq, 
     ylab = "Case Bin Freq (CEU)", 
     xlab = "Control Bin Freq (TSI)",
     xlim=c(0,(max(exon_BF_bins$Control.Freq))), 
     ylim=c(0,(max(top10table$Case.Freq))+.02),
     pch = 20, 
     cex = .3)
title(main = "Case Bin Freq versus Control Bin Freq")
abline(0,1)
text(exon_BF_bins$Control.Freq., exon_BF_bins$Case.Freq., labels_top, cex=0.6, pos=3, col="red")

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