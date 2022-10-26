suppressMessages(library(GGally))
suppressMessages(library(ggtree))
suppressMessages(library(phytools))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Hmisc))
suppressMessages(library(stringr))
suppressMessages(library(reshape2))
suppressMessages(library(reshape))
suppressMessages(library(DESeq2))
suppressMessages(library(gridExtra))
suppressMessages(library(ggforce))
suppressMessages(library(ggrepel))
suppressMessages(library(ggtern))

##########################################################################################
### Virginie Ricci, 2022
##########################################################################################

# Calculate weighted species mean
# Proportion of opsins along the species tree (see Ronco et al. 2021)

##########################################################################################

Method = 'HTSeqCount'
VALUE= 'TPM'

##########################################################################################


##########################################################################################
# RData from HTSeqCount_results.R
##########################################################################################

load(paste0(PATH, Method, '_exons_FilteredForTPMandVST_withRH2As.RData'))


##########################################################################################
# Unpruned species tree (Ronco et al. 2021)
##########################################################################################

species_tree_unpruned_LT_cichlids = read.nexus(paste0(PATH, 'RNAseq_SpeciesTree.tre'))


##########################################################################################
# Total read count per individual (lowly expressed genes removed, protein-coding RNAs and lncRNAs only)
##########################################################################################

ReadCountTotal = data.frame(Species_ID_ID=rownames(HTSeqCount_ALL_ID_init_exons), TotalReadCount=rowSums(HTSeqCount_ALL_ID_init_exons))
rownames(ReadCountTotal) = NULL


##########################################################################################
# Dataframes from RData
##########################################################################################

COUNT_ALL_ID_info_RodCones = get(paste0(Method, '_COUNT_ALL_ID_info_RodCones', '_', VALUE, '_exons'))
COUNT_ALL_ID_info_Cones_SingleDouble = get(paste0(Method, '_COUNT_ALL_ID_info_Cones_SingleDouble', '_', VALUE, '_exons'))

COUNT_ALL_ID_info_Rod_Cones_Single_Double_tmp = COUNT_ALL_ID_info_RodCones[COUNT_ALL_ID_info_RodCones$Opsin == 'RH1',]
COUNT_ALL_ID_info_Rod_Cones_Single_Double = rbind(COUNT_ALL_ID_info_Rod_Cones_Single_Double_tmp, COUNT_ALL_ID_info_Cones_SingleDouble)
COUNT_ALL_ID_info_Rod_Cones_Single_Double = COUNT_ALL_ID_info_Rod_Cones_Single_Double[order(COUNT_ALL_ID_info_Rod_Cones_Single_Double$Species_ID_ID),]


COUNT_ALL_ID = COUNT_ALL_ID_info_Rod_Cones_Single_Double
COUNT_ALL_ID = left_join(COUNT_ALL_ID, ReadCountTotal)

Opsins_final = left_join(data.frame(Opsin=c('SWS1', 'SWS2B', 'SWS2A', 'RH2B', 'RH2Ab', 'RH2Aa', 'LWS', 'RH1')), Opsins)
Opsins_order = Opsins_final$Opsin
Opsins_cols_order = Opsins_final$Colors

# Weighted species mean according to sequencing depth: each opsin out of all opsins
COUNT_ALL_ID_WM = ''
for (sp in unique(COUNT_ALL_ID$Species_ID)){
  totals=c()
  sub_sp = COUNT_ALL_ID[COUNT_ALL_ID$Species_ID== sp,]
  for (ops in unique(COUNT_ALL_ID$Opsin)){
    sub_ops = sub_sp[sub_sp$Opsin == ops,]
    wm = weighted.mean(sub_ops$values_plot, sub_ops$TotalReadCount) # weighted species mean of TPM values corrected for sequencing depth (total read count)
    if (is.na(wm)){
      wm = 0
      print('stop')
    }
    totals = append(totals, wm) # weighted species mean for each opsin
  }
  total = sum(totals) # sum of all weighted species mean for each opsin
  percent = (totals/total)*100
  sum_percent = sum(percent)
  
  TMP = data.frame(Species_ID = '', Opsin = unique(COUNT_ALL_ID$Opsin), GeneName = unique(COUNT_ALL_ID$GeneName), Type='RodCones', WeightedMeanPercent = percent,
                   Tribe = '', Food = '', Habitat = '', Depth = '')
  
  TMP$Species_ID = sp
  TMP$Tribe = unique(sub_ops$Tribe)
  TMP$Food = unique(sub_ops$Food)
  TMP$Habitat = unique(sub_ops$Habitat)
  TMP$Depth = unique(sub_ops$Depth)
  TMP$Tribe_final = unique(sub_ops$Tribe_final)
  
  if (is.data.frame(COUNT_ALL_ID_WM) == FALSE){
    COUNT_ALL_ID_WM = TMP
  }else{
    COUNT_ALL_ID_WM = rbind(COUNT_ALL_ID_WM, TMP)
  }
}
COUNT_ALL_ID_WM_Individually = COUNT_ALL_ID_WM

SINGLECONES = c('SWS1', 'SWS2B', 'SWS2A')
DOUBLCONES = c('RH2B', 'RH2Ab', 'RH2Aa', 'LWS')
COUNT_ALL_ID_WM[COUNT_ALL_ID_WM$Opsin %in% SINGLECONES, ]$Type = 'SingleCones'
COUNT_ALL_ID_WM[COUNT_ALL_ID_WM$Opsin %in% DOUBLCONES, ]$Type = 'DoubleCones'


# Weighted species mean according to sequencing depth: RH1 out of all opsins + SWS1, SWS2B and SWS2A out of single cone opsins + RH2B, RH2Ab, RH2Aa and LWS out of double cone opsins
COUNT_ALL_ID_WM = ''
for (sp in unique(COUNT_ALL_ID$Species_ID)){
  totals=c()
  sub_sp = COUNT_ALL_ID[COUNT_ALL_ID$Species_ID== sp,]
  for (ops in unique(COUNT_ALL_ID$Opsin)){
    sub_ops = sub_sp[sub_sp$Opsin == ops,]
    wm = weighted.mean(sub_ops$values_plot, sub_ops$TotalReadCount) # weighted species mean of TPM values corrected for sequencing depth (total read count)
    if (is.na(wm)){
      wm = 0
      print('stop')
    }
    totals = append(totals, wm) # weighted species mean for each opsin
  }
  names(totals) = unique(COUNT_ALL_ID$Opsin)
  
  total_SingleCones = sum(totals[SINGLECONES]) # sum of all weighted species mean for each opsin in single cones
  total_DoubleCones = sum(totals[DOUBLCONES]) # sum of all weighted species mean for each opsin in double cones
  
  percent_SingleCones = (totals[SINGLECONES]/total_SingleCones)*100
  percent_DoubleCones = (totals[DOUBLCONES]/total_DoubleCones)*100
  
  sum_percent_SingleCones = sum(percent_SingleCones)
  sum_percent_DoubleCones = sum(percent_DoubleCones)
  
  percent = c(percent_SingleCones, percent_DoubleCones)
  percent_order = percent[unique(COUNT_ALL_ID$Opsin)]
  names(percent_order)[1] = unique(COUNT_ALL_ID$Opsin)[1]
  RH1_value = COUNT_ALL_ID_WM_Individually[(COUNT_ALL_ID_WM_Individually$Species_ID == sp) & (COUNT_ALL_ID_WM_Individually$Opsin == 'RH1'),]$WeightedMeanPercent
  percent_order[1] = RH1_value
  
  TMP = data.frame(Species_ID = '', Opsin = names(percent_order), GeneName = unique(COUNT_ALL_ID$GeneName), Type='', WeightedMeanPercent = percent_order,
                   Tribe = '', Food = '', Habitat = '', Depth = '')
  
  TMP$Species_ID = sp
  TMP[TMP$Opsin %in% SINGLECONES,]$Type = 'SingleCones'
  TMP[TMP$Opsin %in% DOUBLCONES,]$Type = 'DoubleCones'
  TMP[TMP$Opsin == 'RH1',]$Type = 'RodCones'
  TMP$Tribe = unique(sub_ops$Tribe)
  TMP$Food = unique(sub_ops$Food)
  TMP$Habitat = unique(sub_ops$Habitat)
  TMP$Depth = unique(sub_ops$Depth)
  TMP$Tribe_final = unique(sub_ops$Tribe_final)
  
  if (is.data.frame(COUNT_ALL_ID_WM) == FALSE){
    COUNT_ALL_ID_WM = TMP
  }else{
    COUNT_ALL_ID_WM = rbind(COUNT_ALL_ID_WM, TMP)
  }
}
COUNT_ALL_ID_WM_Single_Double_Cones = COUNT_ALL_ID_WM







