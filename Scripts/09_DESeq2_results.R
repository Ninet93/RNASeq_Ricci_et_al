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
suppressMessages(library(cluster))
suppressMessages(library(factoextra))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(grid))
suppressMessages(library(EnhancedVolcano))

##########################################################################################
### Virginie Ricci, 2022
##########################################################################################

# Principal Component Analysis (PCA)

##########################################################################################

Method = 'HTSeqCount'

##########################################################################################


##########################################################################################
# Nile tilapia genome annotation
##########################################################################################

FeatureCounts_Orenil_Annot = read.csv(paste0(PATH, 'GCF_001858045.2_O_niloticus_UMD_NMBU_genomic_geneID.gtf_FeatureCount_annotation_exons.txt'), sep='\t', header=T)

Onil_annot_BT = read.csv(paste0(PATH, 'GCF_001858045.2_O_niloticus_UMD_NMBU_genomic_biotypes.txt'), sep='\t', header=T)

Onil_annot_BT_proteincoding_lncRNA = subset(Onil_annot_BT, BioType == 'protein_coding' | BioType == 'lncRNA')


##########################################################################################
# Read count table, TPM, DESeq objects
##########################################################################################

load(paste0(path_RNA, Method, '_exons_FilteredForTPMandVST_withRH2As.RData'))

load(paste0(path_RNA, 'DESeq2/', Method, '_DESeq_output_exons_ProteinCoding_lncRNA_withRH2As.RData'))


##########################################################################################
# Total read count per individual (lowly expressed genes removed, protein-coding RNAs and lncRNAs only)
##########################################################################################

ReadCountTotal = data.frame(Species_ID_ID=rownames(HTSeqCount_ALL_ID_init_exons), TotalReadCount=rowSums(HTSeqCount_ALL_ID_init_exons))
rownames(ReadCountTotal) = NULL



##########################################################################################
# DESeq2 normalized read count table 
##########################################################################################

ddsMethod_Prot_lnc_filter_norm_assay_OPSINS = ddsMethod_Prot_lnc_filter_norm_assay[Opsins$GeneName[-9],]

ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t = data.frame(t(ddsMethod_Prot_lnc_filter_norm_assay_OPSINS))
ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t$Species_ID_ID = rownames(ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t)

ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t_Totals = left_join(ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t, ReadCountTotal)
ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t_Totals$Species_ID = data.frame(do.call(rbind, strsplit(ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t_Totals$Species_ID_ID, '_', 1)))[[1]]


# Weighted species mean on DESeq2 normalized read count table of opsins only
ddsMethod_SpeciesWeightedMean=NULL
for (ops in Opsins$GeneName[c(1:8)]){
  count=count+1
  
  tmp_ddsMethod_SpeciesWeightedMean = as.data.frame(cbind(by(ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t_Totals,
                                                             ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t_Totals$Species_ID,
                                                             function(dat) weighted.mean(dat[,ops], dat$TotalReadCount))))
  names(tmp_ddsMethod_SpeciesWeightedMean) = ops
  
  if (count == 1){
    ddsMethod_SpeciesWeightedMean = tmp_ddsMethod_SpeciesWeightedMean
  }else{
    ddsMethod_SpeciesWeightedMean = cbind(ddsMethod_SpeciesWeightedMean, tmp_ddsMethod_SpeciesWeightedMean)
  }
  
}
ddsMethod_SpeciesWeightedMean$Species_ID = rownames(ddsMethod_SpeciesWeightedMean)


# PCA of weighted species mean values # ROD AND CONE OPSIN GENES
ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t_Totals_split = split(ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t_Totals, ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t_Totals$Species_ID)

Nspecies = length(unique(ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t_Totals$Species_ID))
Nopsins = length(Opsins$GeneName[c(1:8)])
out = as.data.frame(matrix(NA, Nspecies, Nopsins + 1))
out[,1] = unique(ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t_Totals$Species_ID)
names(out)[1] = 'Species_ID'
for (i in 1:Nopsins){
  tmp = lapply(ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t_Totals_split, function(x) weighted.mean(x[,i], x$TotalReadCount))
  out[,i+1] = unlist(tmp)
  names(out)[i+1] = Opsins$GeneName[i]
}
rownames(out) = out$Species_ID
out$Species_ID = NULL

pca_opsins <- prcomp(out)



# PCA of weighted species mean values # CONE OPSIN GENES
ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t_Totals_split = split(ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t_Totals, ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t_Totals$Species_ID)

Nspecies = length(unique(ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t_Totals$Species_ID))
Nopsins = length(Opsins$GeneName[c(2:8)])
out_coneopsins = as.data.frame(matrix(NA, Nspecies, Nopsins + 1))
out_coneopsins[,1] = unique(ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t_Totals$Species_ID)
names(out_coneopsins)[1] = 'Species_ID'
for (i in 1:Nopsins){
  tmp = lapply(ddsMethod_Prot_lnc_filter_norm_assay_OPSINS_t_Totals_split, function(x) weighted.mean(x[,i+1], x$TotalReadCount))
  out_coneopsins[,i+1] = unlist(tmp)
  names(out_coneopsins)[i+1] = Opsins$GeneName[i+1]
}
rownames(out_coneopsins) = out_coneopsins$Species_ID
out_coneopsins$Species_ID = NULL

pca_coneopsins <- prcomp(out_coneopsins)



# Weighted species mean on DESeq2 normalized read count table
ddsMethod_Prot_lnc_filter_norm_assay_t = data.frame(t(ddsMethod_Prot_lnc_filter_norm_assay))
ddsMethod_Prot_lnc_filter_norm_assay_t$Species_ID_ID = rownames(ddsMethod_Prot_lnc_filter_norm_assay_t)

ddsMethod_Prot_lnc_filter_norm_assay_t_Totals = left_join(ddsMethod_Prot_lnc_filter_norm_assay_t, ReadCountTotal)
ddsMethod_Prot_lnc_filter_norm_assay_t_Totals$Species_ID = data.frame(do.call(rbind, strsplit(ddsMethod_Prot_lnc_filter_norm_assay_t_Totals$Species_ID_ID, '_', 1)))[[1]]

ddsMethod_Prot_lnc_filter_norm_assay_t_Totals_split = split(ddsMethod_Prot_lnc_filter_norm_assay_t_Totals, ddsMethod_Prot_lnc_filter_norm_assay_t_Totals$Species_ID)
NcolGenesTotal = nrow(ddsMethod_Prot_lnc_filter_norm_assay)

breaks = seq(1, NcolGenesTotal, 100)
breaks = breaks[-length(breaks)]
breaks_stop = breaks - 1
breaks_stop = breaks_stop[-1]
breaks_stop[length(breaks_stop)+1] = NcolGenesTotal#-1
length(breaks) == length(breaks_stop)
i=1
breaks_i = breaks[i]
breaks_stop_i = breaks_stop[i]

for (i in c(1:length(breaks))){
  print(i)
  breaks_i = breaks[i]
  breaks_stop_i = breaks_stop[i]
  len_breaks = breaks_stop_i-breaks_i+1
  tmp2 = lapply(ddsMethod_Prot_lnc_filter_norm_assay_t_Totals_split, function(x) { apply(x[c(breaks_i:breaks_stop_i)], 2, function(y){ weighted.mean(y, x$Total)})})
  
  out_all2 = as.data.frame(matrix(unlist(tmp2), Nspecies, len_breaks, byrow=T))
  names(out_all2) = colnames(ddsMethod_Prot_lnc_filter_norm_assay_t)[c(breaks_i:breaks_stop_i)]
  if (i == 1){
    out_all = out_all2
  }else{
    out_all = cbind(out_all, out_all2)
  }
  
}

rownames(out_all) = names(ddsMethod_Prot_lnc_filter_norm_assay_t_Totals_split)

pca_all <- prcomp(out_all)


# Protein-coding RNAs
Onil_annot_BT_proteincoding = Onil_annot_BT_proteincoding_lncRNA[Onil_annot_BT_proteincoding_lncRNA$BioType == 'protein_coding',]

PROTEINS = which(colnames(out_all) %in% unique(Onil_annot_BT_proteincoding$GeneID))

ddsMethod_Prot_lnc_filter_norm_assay_PROTEINS = out_all[,PROTEINS]


# lncRNAs ONLY
Onil_annot_BT_lncRNA = Onil_annot_BT_proteincoding_lncRNA[Onil_annot_BT_proteincoding_lncRNA$BioType == 'lncRNA',]

lncRNA = which(colnames(out_all) %in% unique(Onil_annot_BT_lncRNA$GeneID))

ddsMethod_Prot_lnc_filter_norm_assay_lncRNAS = out_all[,lncRNA]




# Phylospace on the PCA with cone opsin genes

phy = read.nexus('Data/RNAseq_SpeciesTree.tre')
# Species tree from Ronco et al. 2021 pruned to taxa included in this study (without Astbur, Oretan, and Tylpol)

# PC1 and PC2
X_df = data.frame('PC1' = pca_coneopsins$x[,1], 'PC2' = pca_coneopsins$x[,2]); head(X_df); dim(X_df)
rownames(X_df) = out_coneopsins$Species_ID; head(X_df); dim(X_df)

points_Tylpol_Oretan_Astbur = X_df[X_df$Species_ID %in% c('Tylpol', 'Oretan', 'Astbur'),]

phylomorphospace(phy, X=X_df, lwd=3, fsize=1.5)
points(points_Tylpol_Oretan_Astbur$PC1, points_Tylpol_Oretan_Astbur$PC2, bg=cols_Tylpol_Oretan_Astbur$Cols, pch=21, cex=1.3)
text(points_Tylpol_Oretan_Astbur$PC1, points_Tylpol_Oretan_Astbur$PC2, labels=points_Tylpol_Oretan_Astbur$Species_ID, cex=0.75)





