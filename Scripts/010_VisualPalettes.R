suppressMessages(library(phytools))

##########################################################################################
### Virginie Ricci, 2022
##########################################################################################

# Ancestral state reconstruction of visual palettes on the species tree

##########################################################################################


##########################################################################################
# input variables
##########################################################################################

SINGLECONES = c('SWS1', 'SWS2B', 'SWS2A')
DOUBLECONES = c('RH2B', 'RH2Aa', 'RH2Ab', 'LWS')

phy # Species tree from Ronco et al. 2021

TPM_WSM # weighted species mean TPM values per species
# species as row and genes as column

TPM_WSM_opsin = TPM_WSM[Opsins[2:8,]$GeneName]  # Opsins table from Orenil_opsins.txt

##########################################################################################
# clustering by gene expression levels
# the single-most expressed cone opsins and the two most expressed double cone opsins
##########################################################################################

SINGLECONES = c('SWS1', 'SWS2B', 'SWS2A')
DOUBLECONES = c('RH2B', 'RH2As', 'LWS')

visual_palettes = data.frame(matrix(ncol=4, nrow=nrow(TPM_WSM_opsin)))
names(visual_palettes) = c('Species', 'Top1', 'Top2', 'Top3')
for (r in seq(1, nrow(TPM_WSM_opsin))){
  TPM_WSM_opsin_r = TPM_WSM_opsin[r,]
  
  id = rownames(TPM_WSM_opsin_r)
  visual_palettes[r,1] = id
  
  t_TPM_WSM_opsin_r = data.frame(t(TPM_WSM_opsin_r))
  rownames(t_TPM_WSM_opsin_r) = Opsins$Opsin[c(2:8)]
  names(t_TPM_WSM_opsin_r) = 'val'
  
  tmp = data.frame('val' = t_TPM_WSM_opsin_r['RH2Aa',] + t_TPM_WSM_opsin_r['RH2Ab',]) # combined RH2Aa and RH2Ab values
  rownames(tmp) = 'RH2As'
  
  t_TPM_WSM_opsin_r = rbind(t_TPM_WSM_opsin_r, tmp)
  
  t_TPM_WSM_opsin_r$opsins = rownames(t_TPM_WSM_opsin_r)
  t_TPM_WSM_opsin_r = t_TPM_WSM_opsin_r[t_TPM_WSM_opsin_r$opsins != 'RH2Aa',] # remove RH2Aa
  t_TPM_WSM_opsin_r = t_TPM_WSM_opsin_r[t_TPM_WSM_opsin_r$opsins != 'RH2Ab',] # remove RH2Ab
  
  t_TPM_WSM_opsin_r_sort = t_TPM_WSM_opsin_r[rev(order(t_TPM_WSM_opsin_r$val)),] # sort df according to values
  order_opsins = rownames(t_TPM_WSM_opsin_r)[rev(order(t_TPM_WSM_opsin_r$val))] # get opsins in the right order

  singlecones=0
  doublecones=0
  for (ops in order_opsins){
    if (ops %in% SINGLECONES){
      singlecones=singlecones+1
      if (singlecones == 1){
        visual_palettes[r,2] = ops # get the most expressed single cone opsin
      }
    }
    if (ops %in% DOUBLECONES){
      doublecones=doublecones+1
      if (doublecones == 1){
        visual_palettes[r,3] = ops # get the first most expressed double cone opsin
      }
      if (doublecones == 2){
        visual_palettes[r,4] = ops # get the second most expressed double cone opsin
      }
    }
  }
}
head(visual_palettes)

df_simmap = paste0(visual_palettes$Top1, '_', visual_palettes$Top2, '_', visual_palettes$Top3)
names(df_simmap) = visual_palettes$Species	

df_simmap[df_simmap == 'SWS1_LWS_RH2As'] = 'SWS1_RH2As_LWS'
df_simmap[df_simmap == 'SWS1_LWS_RH2B'] = 'SWS1_RH2B_LWS'
df_simmap[df_simmap == 'SWS1_RH2As_RH2B'] = 'SWS1_RH2B_RH2As'
df_simmap[df_simmap == 'SWS2A_LWS_RH2As'] = 'SWS2A_RH2As_LWS'
df_simmap[df_simmap == 'SWS2A_RH2As_RH2B'] = 'SWS2A_RH2B_RH2As'
df_simmap[df_simmap == 'SWS2B_LWS_RH2As'] = 'SWS2B_RH2As_LWS'
df_simmap[df_simmap == 'SWS2B_LWS_RH2B'] = 'SWS2B_RH2B_LWS'
df_simmap[df_simmap == 'SWS2B_RH2As_RH2B'] = 'SWS2B_RH2B_RH2As'
levels(df_simmap) = c('SWS1_RH2As_LWS','SWS1_RH2B_LWS','SWS1_RH2B_RH2As',
  'SWS2A_RH2As_LWS','SWS2A_RH2B_RH2As',
  'SWS2B_RH2As_LWS', 'SWS2B_RH2B_LWS','SWS2B_RH2B_RH2As')


##########################################################################################
# clustering by hierarchical clustering method (Ward's method)
##########################################################################################

TPM_WSM_DoubleSingleCones = TPM_WSM_opsin

TPM_WSM_DoubleSingleCones[SINGLECONES] = TPM_WSM_DoubleSingleCones[SINGLECONES] / rowSums(TPM_WSM_DoubleSingleCones[SINGLECONES])
TPM_WSM_DoubleSingleCones[DOUBLCONES] = TPM_WSM_DoubleSingleCones[DOUBLCONES] / rowSums(TPM_WSM_DoubleSingleCones[DOUBLCONES])

hc = hclust(dist(TPM_WSM_DoubleSingleCones), method = 'ward.D2')

df_simmap_hclust = cutree(hc, k = 3)
df_simmap_hclust[df_simmap_hclust == 1] = 'SWS2B_RH2As_RH2B'
df_simmap_hclust[df_simmap_hclust == 2] = 'SWS1_RH2As_RH2B'
df_simmap_hclust[df_simmap_hclust == 3] = 'SWS2A_LWS_RH2As'
levels(df_simmap_hclust) = c('SWS1_RH2As_RH2B', 'SWS2B_RH2As_RH2B', 'SWS2A_LWS_RH2As')

df_simmap_hclust[df_simmap_hclust == 'SWS1_RH2As_RH2B'] = 'SWS1_RH2B_RH2As'
df_simmap_hclust[df_simmap_hclust == 'SWS2A_LWS_RH2As'] = 'SWS2A_RH2As_LWS'
df_simmap_hclust[df_simmap_hclust == 'SWS2B_RH2As_RH2B'] = 'SWS2B_RH2B_RH2As'
levels(df_simmap_hclust) = c('SWS1_RH2B_RH2As', 'SWS2A_RH2As_LWS', 'SWS2B_RH2B_RH2As')


##########################################################################################
# ancestral state reconstruction with make.simmap
##########################################################################################

parameters=c(p_type, p_model, p_method, p_nsim)
names(parameters) = c('type', 'model', 'method', 'nsim')

mtree = make.simmap(phy, df_simmap, type = parameters[['type']], model = parameters[['model']], method = parameters[['method']], nsim = as.numeric(parameters[['nsim']]), Q='mcmc')

mtree_hclust = make.simmap(phy, df_simmap_hclust, type = parameters[['type']], model = parameters[['model']], method = parameters[['method']], nsim = as.numeric(parameters[['nsim']]), Q='mcmc')











