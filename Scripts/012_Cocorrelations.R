suppressMessages(library(phytools))
suppressMessages(library(caper))

##########################################################################################
### Fabrizia Ronco and Virginie Ricci, 2022
##########################################################################################

# Co-correlation of visual opsin gene expression (lm and PIC analyses)

##########################################################################################


##########################################################################################
# input variables
##########################################################################################

phy = read.nexus('Data/RNAseq_SpeciesTree.tre')
# Species tree from Ronco et al. 2021 pruned to taxa included in this study (without Astbur, Oretan, and Tylpol)

Opsins = read.csv('Data/Orenil_opsins', sep='\t', header=F)
names(Opsins) = c('Opsin', 'GeneName', 'Color')

TPM_WSM # output of the script 08_HTSeqCount_results_WeightedSpeciesMean_opsins.R
# == COUNT_ALL_ID_WM_Individually
# weighted species mean TPM values per species (without Astbur, Oretan, and Tylpol)
# species as row and genes as column


##########################################################################################
# Expression correlation of visual opsin gene expression
##########################################################################################

TPM_WSM_opsin = TPM_WSM[Opsins$Opsin,]

TPM_WSM_opsin$RH2As = TPM_WSM_opsin$RH2Aa + TPM_WSM_opsin$RH2Ab

TPM_WSM_opsin_Total = apply(TPM_WSM_opsin, 1, sum)

TPM_WSM_opsin_Percent = TPM_WSM_opsin

TPM_WSM_opsin_Percent = TPM_WSM_opsin_Percent / TPM_WSM_opsin_Total
# colnames = RH1, SWS1, SWS2A, SWS2B, RH2B, RH2Aa, RH2Ab, LWS, RH2As

PICS = matrix(NA, nrow(TPM_WSM_opsin_Percent)-1, ncol(TPM_WSM_opsin_Percent))
colnames(PICS)=colnames(TPM_WSM_opsin_Percent)
for (i in 1:ncol(TPM_WSM_opsin_Percent)) {
	ops = colnames(TPM_WSM_opsin_Percent)[i]
	tmp = TPM_WSM_opsin_Percent[,ops]; names(tmp)=TPM_WSM_opsin_Percent$Species_ID
	PICS[,i] = pic(tmp, phy)
}

cor_res = cor.table(PICS, cor.method="spearman", cor.type="contrast")
cor_res_raw = cor.table(TPM_WSM_opsin_Percent, cor.method="spearman", cor.type="standard")


par(mfrow=c(9,9), mar=c(0.5,0.5,0.5,0.5))
for (i in 1:ncol(TPM_WSM_opsin_Percent)) {
  genei = colnames(TPM_WSM_opsin_Percent)[i]
  for (k in 1:ncol(TPM_WSM_opsin_Percent)) {
    genek = colnames(TPM_WSM_opsin_Percent)[k]
    if (i == k) {
      plot(0,0,type="n", axes=F, xlab="", ylab="")
      text(0,0, colnames(PICS)[k], cex=2)
    }
    if (i < k) {
      plot(TPM_WSM_opsin_Percent[,genei]~TPM_WSM_opsin_Percent[,genek], las=1, axes=F, pch=19, cex=0.75)
      box() 
      if ( k-i == 1 ) {
        axis(2, tck=-0.02, las=1, cex.axis=0.5, hadj=0.3)
        axis(1, tck=-0.02, las=1, cex.axis=0.5, padj=-4)
      } else {
        axis(2, tck=-0.02, labels=F)
        axis(1, tck=-0.02, labels=F)
      }
      mod = lm(TPM_WSM_opsin_Percent[,genei] ~ TPM_WSM_opsin_Percent[,genek])
      legend("topright", paste("r =",round(cor_res_raw$r[k,i],2), ";p = ", round(cor_res_raw$P[k,i],4)) ,bty="n", cex=0.75)	
      lty="solid"
      if (cor_res_raw$P[k,i] >= 0.05 ) { lty="dashed" }
      abline(mod, col="grey40", lty=lty)
    }
    if (i > k) { 
      plot(PICS[,i]~PICS[,k], col="grey30", las=1, axes=F, pch=19, cex=0.75)
      box() 
      if ( i-k == 1 ) {
        axis(4, tck=-0.02, las=1, cex.axis=0.5, hadj=0.8)
        axis(3, tck=-0.02, las=1, cex.axis=0.5, padj=3)
      } else {
        axis(4, tck=-0.02, labels=F)
        axis(3, tck=-0.02, labels=F)
      }
      mod = lm(PICS[,i] ~ PICS[,k] - 1 )
      legend("bottomleft", paste("r =",round(cor_res$r[k,i],2), ";p = ", round(cor_res$P[k,i],4)) ,bty="n", cex=0.75)	
      lty="solid"
      if ( cor_res$P[k,i] >= 0.05 ) { lty="dashed" }
      abline(mod, col="grey40", lty=lty)
    }
    
  }
}



### P-values multiple testing correction
# including both RH2Aa, RH2Ab and RH2As

Ps = cor_res$P
diag(Ps)=NA
Ps[upper.tri(Ps)]=NA
Ps = melt(Ps, na.rm=T)
Ps$corrected = p.adjust(Ps$value, method="BH")
write.csv(Ps, "P_values_corrected_9genes_Spearman_PICs.csv", row.names=F)


Ps_raw=cor_res_raw$P
diag(Ps_raw)=NA
Ps_raw[upper.tri(Ps_raw)]=NA
Ps_raw = melt(Ps_raw, na.rm=T)
Ps_raw$corrected = p.adjust(Ps_raw$value, method="BH")
write.csv(Ps_raw, "P_values_corrected_9genes_Spearman_raw.csv", row.names=F)

