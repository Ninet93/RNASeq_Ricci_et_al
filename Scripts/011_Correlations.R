suppressMessages(library(phytools))
suppressMessages(library(geomorph))

##########################################################################################
### Virginie Ricci, 2022
##########################################################################################

# Linear regression (lm) and phylogenetic generalized least squares (pGLS) 

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

out # output of the script 09_DESeq2_results.R
# dds of opsins
# VST normalized values per species
# species as row and opsin genes as column


##########################################################################################
# lm and pGLS between opsin expression and eco-morphological traits
##########################################################################################

Traits = c('d13N', 'd15N', 'EyeSize')

### lm and pGLS results
stats_df = data.frame(matrix(ncol=7, nrow=0))
names(stats_df) = c('Trait', 'Opsin', 'lm_R2', 'lm_pval', 'pGLS_lambda', 'pGLS_R2', 'pGLS_pval')
###

for (TRAIT in Traits){
for (OPSIN in Opsins$Opsin){
  TPM_WSM_opsin = TPM_WSM[TPM_WSM$Opsin == OPSIN,]
  TPM_WSM_opsin_trait = subset(TPM_WSM_opsin, select=c('Species_ID', 'WeightedMeanVST', TRAIT))
  
  ### lm and pGLS summary
  formula_df = TPM_WSM_opsin_trait[, 'WeightedMeanVST'] ~ TPM_WSM_opsin_trait[, TRAIT]
  lm_summary = summary(lm(formula_df))
  lm_rsquared = round(lm_summary$r.squared, 4)
  lm_pval = as.numeric(round(lm_summary$coefficients[,4][2], 4))
  
  compset = comparative.data(phy, TPM_WSM_opsin_trait, Species_ID)
  pgls_summary = summary(pgls(WeightedMeanVST ~ TRAIT, compset, lambda='ML'))
  pgls_lambda = as.numeric(round(pgls_summary$param[2], 4))
  pgls_rsquared = round(pgls_summary$r.squared, 4)
  pgls_pval = as.numeric(round(pgls_summary$coefficients[,4][2], 4))
  ### lm and pGLS summary
  
  stats_df_toadd = data.frame(X=TRAIT, Opsin=OPSIN, lm_R2=lm_rsquared, lm_pval=lm_pval, pGLS_lambda=pgls_lambda, pGLS_R2=pgls_rsquared, pGLS_pval=pgls_pval)
  stats_df = rbind(stats_df,stats_df_toadd)
  
  }
}


##########################################################################################
# lm and two-block PLS between cone opsin expression and eco-morphological traits
##########################################################################################

VST_WSM_coneopsins_trait # merge of out and table 1 of the article, select cone opsins only
# VST values per species with d13C, d15N and eye size values
# species as row

for (TRAIT in Traits){
  plsout = phylo.integration(VST_WSM_coneopsins_trait[TRAIT], VST_WSM_coneopsins_trait[cone_opsins], phy)
  
  f1= data.frame('XScores'= plsout$XScores[,1], 'YScores'= plsout$YScores[,1])
  mod = summary(lm(f1$YScores ~ f1$XScores))
  
  plot(f1$XScores, f1$YScores, pch=19, las=1, xlab=TRAIT, ylab='Cone opsins')
  abline(lm(f1$YScores ~ f1$XScores))
  legend( 'topleft', paste0('r.pls=',round(plsout$r.pls,3) ,', p=', round(plsout $P.value,3), '\nR2=',round(mod$r.squared,3) ,', p=', round(mod$coefficients[2,4],3) ) )
  legend('bottomright', paste0( rownames(plsout$right.pls.vectors), '     ' , round(plsout$right.pls.vectors, 2)), cex=0.75)
  legend('bottomleft', paste0( rownames(plsout$left.pls.vectors), '     ' , round(plsout$left.pls.vectors, 2)), cex=0.75)

}