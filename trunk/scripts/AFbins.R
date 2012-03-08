exon_AF_bins <- read.csv("/Users/carriecb/Desktop/BioBin/TSI.CEUcomp/Add_comp/MAF.05_Add_gene-only_tests/combined_output_b37-AllFreq.csv")
#to get rid of common variants
exon_AF_rare_bins <- exon_AF_bins[0,]
for (i in 1:length(exon_AF_bins$Locus)){
  if (exon_AF_bins[i,4] == 1){
    exon_AF_rare_bins <- rbind(exon_AF_rare_bins, exon_AF_bins[i,])
  }
}

resid <- abs(exon_AF_rare_bins$Control.MAF - exon_AF_rare_bins$Case.MAF)/sqrt(2)
exon_AF_rare_bins <- cbind(exon_AF_rare_bins, resid)
for (i in 1:length(exon_AF_rare_bins$Locus)){
  if (exon_AF_rare_bins[i,2]<0 || exon_AF_rare_bins[i,3]<0){
    exon_AF_rare_bins[i,6] <- -1
  }

}
View(exon_AF_rare_bins)
#find the limit of the highest ten residuals
sorted_resid <- sort(exon_AF_rare_bins$resid, TRUE)
resid_lim <- sorted_resid[30]

#create a column with ids for only those top 10 residuals
labels_top <- exon_AF_rare_bins$Locus
top10table <- exon_AF_rare_bins[0,]
for (i in 1:length(exon_AF_rare_bins$Locus)){
  if (exon_AF_rare_bins[i,6] < resid_lim){
    labels_top[i] <- ''
  }
  else{
    top10table <- rbind(top10table,exon_AF_rare_bins[i,])
  }
}
show(top10table)
# top10table1 <- as.data.frame(top10table)
# View(top10table)
# sortedtop10table <- top10table1[order(resid) , ]
# View(sortedtop10table)

plot(exon_AF_rare_bins$Control.MAF, exon_AF_rare_bins$Case.MAF,
     ylab = "Case MAF (CEU)",
     xlab = "Control MAF (TSI)",
     xlim=c(0,.05),
     ylim=c(0,(max(top10table$Case.MAF))+.02),
     pch = 20,
     cex = .3)
title(main = "Case MAF versus Control MAF")
abline(0,1)
text(exon_AF_rare_bins$Control.MAF, exon_AF_rare_bins$Case.MAF, labels_top, cex=0.6, pos=3, col="blue")
