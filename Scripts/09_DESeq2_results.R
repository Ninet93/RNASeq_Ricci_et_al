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


# PCA of weighted species mean values
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






