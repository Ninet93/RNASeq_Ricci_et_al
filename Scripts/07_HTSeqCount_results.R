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
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

##########################################################################################
### Virginie Ricci, 2022
##########################################################################################

# Create read count table from HTSeq-count outputs
# Select expressed protein-coding RNAs and lncRNAs 
# Remove lowly expressed genes
# Convert read count to TPM
# Convert TPM to %TPM per opsin per opsin category
# Create DESeq2 object

##########################################################################################

Method = 'HTSeqCount'
DESeq_TF = TRUE

##########################################################################################


##########################################################################################
# TPM normalisation - counts per length of transcript (kb) per million reads mapped 
# account for sequencing depth and gene length
##########################################################################################

tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}


##########################################################################################
# Stats per individual
##########################################################################################

Stats_per_row = function(x_row){
  ID = rownames(x_row)
  total_row=rowSums(x_row)
  ExpressedGenes = length(which(x_row != 0))
  ExpressedGenes_Percent = (ExpressedGenes/ncol(x_row))*100
  df = data.frame(Species_ID_ID=ID, TotalCount_TPM=total_row, ExpressedGenes_TPM=ExpressedGenes, ExpressedGenes_TPM_Percent=ExpressedGenes_Percent)
  return(df)
}

Stats_df = function(x){
  min_t = rowQuantiles(as.matrix(x))[,1]
  quart1_t = rowQuantiles(as.matrix(x))[,2]
  median_t = rowQuantiles(as.matrix(x))[,3]
  mean_t = rowMeans(as.matrix(x))
  quart3_t = rowQuantiles(as.matrix(x))[,4]
  max_t = rowQuantiles(as.matrix(x))[,5]
  
  df = data.frame(MinCountInGenes_TPM=min_t, Quart1CountInGenes_TPM=quart1_t, MedianCountInGenes_TPM=median_t,
                  MeanCountInGenes_TPM=mean_t, Quart3CountInGenes_TPM=quart3_t, MaxCountInGenes_TPM=max_t)
  return(df)
}


##########################################################################################
# Opsin information
##########################################################################################

Opsins = read.csv(paste0(PATH, 'Orenil_opsins.txt'), sep='\t', header=FALSE)
names(Opsins) = c('Opsin', 'GeneName', 'Colors')


##########################################################################################
# Nile tilapia genome annotation
##########################################################################################

FeatureCounts_Orenil_Annot = read.csv(paste0(PATH, 'GCF_001858045.2_O_niloticus_UMD_NMBU_genomic_geneID.gtf_FeatureCount_annotation_exons.txt'), sep='\t', header=T)

Onil_annot_BT = read.csv(paste0(PATH, 'GCF_001858045.2_O_niloticus_UMD_NMBU_genomic_biotypes.txt'), sep='\t', header=T)

Onil_annot_BT_proteincoding_lncRNA = subset(Onil_annot_BT, BioType == 'protein_coding' | BioType == 'lncRNA')


##########################################################################################
# Cichlid species information, colors for figures, summary statistics
##########################################################################################

# Tribe colors
Tribe_cols = unique(subset(Df_plot_final, select=c(Tribe, Cols)))


# Species information (from Ricci et al. 2022)
Infos = read.csv(paste0(PATH, 'RH1_individuals_info_article_food_habitat_depth.txt'), sep='\t', header=TRUE)
Infos$Species_ID = tolower(Infos$Species_ID)
substring(Infos$Species_ID, 1, 1) = toupper(substring(Infos$Species_ID, 1, 1))
Infos = rbind(Infos, data.frame(Species_ID='Tylpol', Tribe='Tylochromini', Food=NA, Habitat=NA, Depth=NA))

Infos_cols = Infos
Infos_cols[is.na(Infos_cols)] = ''
Infos_cols[Infos_cols$Habitat == 'inter',]$Habitat = 'intermediate'


# Habitat colors
colfunc<-colorRampPalette(c("gold", "limegreen", "olivedrab4", "royalblue3", "mediumorchid4"))
col_habitat = colfunc(length(unique(Infos_cols$Habitat)))

Infos_cols$Habitat_colors = ''
count=0
for (var in c('littoral', 'shallow', 'pelagic', 'intermediate', 'deep')){
  count = count + 1
  Infos_cols[Infos_cols$Habitat == var,]$Habitat_colors = col_habitat[count]
}


# Food colors
colfunc<-colorRampPalette(c("red", "gold", "limegreen", "olivedrab4", "lightskyblue", "royalblue3", "magenta4", "hotpink"))
col_food = colfunc(length(unique(Infos_cols$Food)))

Infos_cols$Food_colors = ''
count=0
for (var in c('fish', 'fish_invert', 'invert', 'auf_invert', 'auf_herb', 'herb', 'plan_invert', 'plan', 'fry_plan', 'scales', 'omni')){
  count = count + 1
  Infos_cols[Infos_cols$Food == var,]$Food_colors = col_food[count]
}


# Depth colors
Infos_cols$Depth_colors = ''
for (var in c('intermediate', 'shallow', 'deep')){
  Infos_cols[Infos_cols$Depth == var,]$Depth_colors = unique(Infos_cols[Infos_cols$Habitat == var,]$Habitat_colors)
}


# Sampling information
Infos_LabKey = read.csv(paste0(PATH, SAMPLING_INFO), header=TRUE)
names(Infos_LabKey) = c('Species_ID', 'Species', 'ID', 'TissueTubeID', 'Sex', 'Tribe', 'Origin', 'Location', 'Date', 'TissueID', 'Notes', 'Morph')
Infos_LabKey$Species_ID = tolower(Infos_LabKey$Species_ID)
substring(Infos_LabKey$Species_ID, 1, 1) = toupper(substring(Infos_LabKey$Species_ID, 1, 1))
Infos_LabKey$Species_ID_ID = paste0(Infos_LabKey$Species_ID, '_', str_replace(Infos_LabKey$ID, 'S_', ''), '_', str_replace(Infos_LabKey$TissueTubeID, '-eye', ''))
sub_Infos_LabKey = subset(Infos_LabKey, select=c(Species_ID, Species_ID_ID, Sex, Origin, Location, Date))

Infos_LabKey_final_tmp = left_join(Infos_cols, sub_Infos_LabKey)


# Table with sample IDs and HTSeq-count output file names
head(Method_count)


##########################################################################################
# HTSeq-count output concatenation
##########################################################################################

start_time = Sys.time()
ALL_ID_init0 = ''
# Create dataframe with counts for all individuals
Expression_all_ID = ''
for (ID in unique(Method_count$Species_ID_ID)){
  path_ID=unique(Method_count[grep(ID, Method_count$Species_ID_ID),]$FilePath)
  
  count_ID = as.data.frame(fread(unique(paste0(path_ID, Method_count[grep(ID, Method_count$Species_ID_ID),]$FileName)), sep='\t', header=FALSE))
  names(count_ID) = c('GeneName', 'Count')
  
  count_ID$Species_ID_ID = ID
  count_ID$Species_ID = data.frame(do.call(rbind, strsplit(count_ID$Species_ID_ID, '_', 1)))[,1]
  
  ###
  path_ID_RH2As = unique(Method_count[grep(ID, Method_count$Species_ID_ID),]$FilePathRH2As)
  
  count_ID_RH2As = as.data.frame(fread(unique(paste0(path_ID_RH2As, Method_count[grep(ID, Method_count$Species_ID_ID),]$FileNameRH2As)), sep='\t', header=FALSE))
  names(count_ID_RH2As) = c('GeneName', 'Count')
  
  count_ID_RH2As$Species_ID_ID = ID
  count_ID_RH2As$Species_ID = data.frame(do.call(rbind, strsplit(count_ID_RH2As$Species_ID_ID, '_', 1)))[,1]
  
  # keeping the proportion of RH2Aa and RH2Ab when adding unmapped read counts
  tot_RH2As = count_ID[count_ID$GeneName == 'LOC100710942',]$Count + count_ID[count_ID$GeneName == 'LOC100710676',]$Count
    
  prop_RH2Aa = count_ID[count_ID$GeneName == 'LOC100710942',]$Count/tot_RH2As
  prop_RH2Ab = count_ID[count_ID$GeneName == 'LOC100710676',]$Count/tot_RH2As
    
  tot_unmapped_RH2As = count_ID_RH2As[count_ID_RH2As$GeneName == 'LOC100710942',]$Count + count_ID_RH2As[count_ID_RH2As$GeneName == 'LOC100710676',]$Count
    
  extra_RH2Aa = tot_unmapped_RH2As*prop_RH2Aa
  extra_RH2Ab = tot_unmapped_RH2As*prop_RH2Ab
    
  count_ID_RH2As[count_ID_RH2As$GeneName == 'LOC100710942', ]$Count = extra_RH2Aa
  count_ID_RH2As[count_ID_RH2As$GeneName == 'LOC100710676', ]$Count = extra_RH2Ab
  
  
  ###
  matrixx = as.data.frame(t(count_ID['Count']))
  colnames(matrixx) = count_ID$GeneName
  
  matrixx_RH2As = as.data.frame(t(count_ID_RH2As['Count']))
  colnames(matrixx_RH2As) = count_ID_RH2As$GeneName
  
  ###
  matrixx_complete = matrixx + matrixx_RH2As
  matrixx_complete$Species_ID_ID = ID
  
  matrixx_final = matrixx_complete[c(length(count_ID$GeneName)+1,1:length(count_ID$GeneName))]
  
  matrixx_final_copy = matrixx_final[2:ncol(matrixx_final)]
  total_row = rowSums(as.matrix(matrixx_final_copy))
  ExpressedGenes = length(which(matrixx_final_copy != 0))
  ExpressedGenes_Percent = (ExpressedGenes/length(count_ID$GeneName))*100
  min_t = min(as.matrix(matrixx_final_copy))
  quart1_t = quantile(as.matrix(matrixx_final_copy), 0.25)[[1]]
  median_t = median(as.matrix(matrixx_final_copy))
  mean_t = mean(as.matrix(matrixx_final_copy))
  quart3_t = quantile(as.matrix(matrixx_final_copy), 0.75)[[1]]
  max_t = max(as.matrix(matrixx_final_copy))
  
  Expression_df = data.frame(Species_ID_ID=ID, TotalCount=total_row, ExpressedGenes=ExpressedGenes, ExpressedGenes_Percent=ExpressedGenes_Percent, MinCountInGenes=min_t, Quart1CountInGenes=quart1_t, MedianCountInGenes=median_t, MeanCountInGenes=mean_t, Quart3CountInGenes=quart3_t, MaxCountInGenes=max_t)
  
  if (is.data.frame(ALL_ID_init0) == FALSE){
    ALL_ID_init0 = matrixx_final
    Expression_all_ID = Expression_df
  }else{
    ALL_ID_init0 = rbind(ALL_ID_init0, matrixx_final)
    Expression_all_ID = rbind(Expression_all_ID, Expression_df)
  }
  print(nrow(ALL_ID_init0))
  
}
stop_time = Sys.time()
start_time; stop_time
rm(count_ID, matrixx, matrixx_final)
# I round() values of ALL_ID_init0, otherwise dds doesn't work
ALL_ID_init = round(ALL_ID_init0) # BackUp
rownames(ALL_ID_init) = ALL_ID_init$Species_ID_ID
ALL_ID_init$Species_ID_ID = NULL

print('HTSeq-count output concatenation done...')


##########################################################################################
# Select expressed protein-coding RNAs and lncRNAs and remove lowly expressed genes + DESeq2 objects
##########################################################################################

ALL_ID_filter = ALL_ID_init[which(colnames(ALL_ID_init) %in% Onil_annot_BT_proteincoding_lncRNA$GeneID)]

ddsMethod_df_tmp = left_join(Method_count, Infos_LabKey_final)
ddsMethod_df = unique(ddsMethod_df_tmp[c('Species_ID_ID', 'FileName', 'Species_ID', 'Sex', 'Tribe', 'Food', 'Habitat', 'Depth')])
ddsMethod_df$FolderFileName = paste0(ddsMethod_df$Species_ID_ID, '_trimmed_singlemapping_exons/', ddsMethod_df$FileName)
ddsMethod_df = ddsMethod_df[c('Species_ID_ID', 'FolderFileName', 'Species_ID', 'Sex', 'Tribe', 'Food', 'Habitat', 'Depth')]

dds <- DESeqDataSetFromMatrix(countData = t(ALL_ID_filter),
                              colData = ddsMethod_df,
                              design= ~ Species_ID + Sex)

ALL_ID_filter = t(assay(dds[rowSums(counts(dds)>5) >= 3,]))


ddsMethod_Prot_lnc_filter = ALL_ID_filter

ddsMethod_Prot_lnc_filter_norm = vst(ddsMethod_Prot_lnc_filter, blind=FALSE)
ddsMethod_Prot_lnc_filter_norm_assay = assay(ddsMethod_Prot_lnc_filter_norm)
ddsMethod_Prot_lnc_filter_norm_assay_dist = as.matrix(dist(t(ddsMethod_Prot_lnc_filter_norm_assay)))

dds_Wald_Prot_lnc = DESeq(ddsMethod_Prot_lnc_filter, test="Wald")

save(ddsMethod_Prot_lnc_filter_norm, ddsMethod_Prot_lnc_filter_norm_assay, ddsMethod_Prot_lnc_filter_norm_assay_dist, dds_Wald_Prot_lnc, file=paste0(PATH, Method, '_DESeq_output_exons_ProteinCoding_lncRNA_withRH2As.RData'))


##########################################################################################
# TPM normalisation
##########################################################################################

ALL_ID_filter_T = t(ALL_ID_filter) # samples (columns) and genes (rows)
ALL_ID_TPM = data.frame(t(tpm3(ALL_ID_filter_T, FeatureCounts_Orenil_Annot[FeatureCounts_Orenil_Annot$GeneID %in% colnames(ALL_ID_filter),]$Length)))

print('TPM normalisation done...')


##########################################################################################
# Read count and TPM for opsins only
##########################################################################################

start_time = Sys.time()
# Create dataframe with counts for all individuals - only opsins >> "RH1" "RH2Aa" "RH2Ab" "RH2B" "SWS1" "SWS2A" "SWS2B" "LWS"
opsins_all_ID = data.frame(matrix(ncol=6, nrow=0))
opsins_all_ID_TPM = data.frame(matrix(ncol=6, nrow=0))
for (ID in unique(Method_count$Species_ID_ID)){
  path_ID=unique(Method_count[grep(ID, Method_count$Species_ID_ID),]$FilePath)
  
  count_ID = as.data.frame(fread(unique(paste0(path_ID, Method_count[grep(ID, Method_count$Species_ID_ID),]$FileName)), sep='\t', header=FALSE))
  names(count_ID) = c('GeneName', 'Count')
  count_ID$Species_ID_ID = ID
  count_ID$Species_ID = data.frame(do.call(rbind, strsplit(count_ID$Species_ID_ID, '_', 1)))[,1]
  
  opsins_ID = count_ID[which(count_ID$GeneName %in% Opsins$GeneName),]
  
  opsins_ID = left_join(opsins_ID, Opsins)
  opsins_ID$Type = 'Cones'
  opsins_ID[opsins_ID$Opsin == 'RH1', ]$Type = 'Rod'
  
  opsins_ID$Total = 0
  opsins_ID$Percent = 0
  
  opsins_ID[opsins_ID$Type == 'Cones',]$Total = sum(opsins_ID[opsins_ID$Type == 'Cones',]$Count)
  opsins_ID[opsins_ID$Type == 'Rod',]$Total = sum(opsins_ID[opsins_ID$Type == 'Rod',]$Count)
  opsins_ID$Percent = (opsins_ID$Count / opsins_ID$Total)*100
  
  opsins_ID_final = opsins_ID[c('Species_ID_ID', 'Species_ID', 'Opsin', 'GeneName', 'Type', 'Count', 'Total', 'Percent', 'Colors')]
  opsins_all_ID = rbind(opsins_all_ID, opsins_ID_final)
  print(length(unique(opsins_all_ID$Species_ID_ID)))
  
  ### TPM
  ALL_ID_TPM_ID = data.frame(t(ALL_ID_TPM[ID,]))
  ALL_ID_TPM_ID$GeneName = rownames(ALL_ID_TPM_ID)
  rownames(ALL_ID_TPM_ID) = NULL
  names(ALL_ID_TPM_ID)[1] = 'TPM'
  ALL_ID_TPM_ID$Species_ID_ID = ID
  ALL_ID_TPM_ID$Species_ID = data.frame(do.call(rbind, strsplit(ALL_ID_TPM_ID$Species_ID_ID, '_', 1)))[,1]
  ALL_ID_TPM_ID = ALL_ID_TPM_ID[c('GeneName', 'TPM', 'Species_ID_ID', 'Species_ID')]
  opsins_ID_TPM = ALL_ID_TPM_ID[which(ALL_ID_TPM_ID$GeneName %in% Opsins$GeneName),]
  
  opsins_ID_TPM = left_join(opsins_ID_TPM, Opsins)
  opsins_ID_TPM$Type = 'Cones'
  opsins_ID_TPM[opsins_ID_TPM$Opsin == 'RH1', ]$Type = 'Rod'
  
  opsins_ID_TPM$Total = 0
  opsins_ID_TPM$Percent = 0
  
  opsins_ID_TPM[opsins_ID_TPM$Type == 'Cones',]$Total = sum(opsins_ID_TPM[opsins_ID_TPM$Type == 'Cones',]$TPM)
  opsins_ID_TPM[opsins_ID_TPM$Type == 'Rod',]$Total = sum(opsins_ID_TPM[opsins_ID_TPM$Type == 'Rod',]$TPM)
  opsins_ID_TPM$Percent = (opsins_ID_TPM$TPM / opsins_ID_TPM$Total)*100
  
  opsins_ID_TPM_final = opsins_ID_TPM[c('Species_ID_ID', 'Species_ID', 'Opsin', 'GeneName', 'Type', 'TPM', 'Total', 'Percent', 'Colors')]
  opsins_all_ID_TPM = rbind(opsins_all_ID_TPM, opsins_ID_TPM_final)
  print(length(unique(opsins_all_ID_TPM$Species_ID_ID)))
  
}
stop_time = Sys.time()
start_time; stop_time
rm(count_ID, opsins_ID, opsins_ID_final, opsins_ID_TPM, opsins_ID_TPM_final)

print('Read count and TPM for opsins only done...')


##########################################################################################
# Opsins per category
##########################################################################################

for (COUNT_ALL_ID_n in c('opsins_all_ID', 'opsins_all_ID_TPM')){
  if (COUNT_ALL_ID_n == 'opsins_all_ID'){
    value = 'Count'
    COUNT_ALL_ID = get(COUNT_ALL_ID_n)
    names(COUNT_ALL_ID)[which(names(COUNT_ALL_ID) == 'Count')] = 'values_plot'
  }
  if (COUNT_ALL_ID_n == 'opsins_all_ID_TPM'){
    value = 'TPM'
    COUNT_ALL_ID = get(COUNT_ALL_ID_n)
    names(COUNT_ALL_ID)[which(names(COUNT_ALL_ID) == 'TPM')] = 'values_plot'
  }
}

###
COUNT_ALL_ID_info = unique(left_join(COUNT_ALL_ID, Infos_LabKey_final))
COUNT_ALL_ID_info$Sex_Species_ID_ID = paste0(COUNT_ALL_ID_info$Sex, '_', COUNT_ALL_ID_info$Species_ID_ID)


##########################################################################################
# Cone opsins only
##########################################################################################

COUNT_ALL_ID_info_Cones = COUNT_ALL_ID_info[COUNT_ALL_ID_info$Type == 'Cones', ]

Opsins_tmp = unique(subset(COUNT_ALL_ID_info_Cones, select=c(Opsin,Colors)))
Opsins_final = left_join(data.frame(Opsin=c('SWS1', 'SWS2B', 'SWS2A', 'RH2B', 'RH2Ab', 'RH2Aa', 'LWS')), Opsins_tmp)
Opsins_order = Opsins_final$Opsin
Opsins_cols_order = Opsins_final$Colors


##########################################################################################
# Rhodopsin and cone opsins together
##########################################################################################

COUNT_ALL_ID_info_RodCones = subset(COUNT_ALL_ID_info, Type == 'Cones' | Type == 'Rod')
COUNT_ALL_ID_info_RodCones[COUNT_ALL_ID_info_RodCones$Type != 'exoRod',]$Type = 'RodCones'

for (ID in Method_count$Species_ID_ID){
  COUNT_ALL_ID_info_RodCones[(COUNT_ALL_ID_info_RodCones$Species_ID_ID == ID) & (COUNT_ALL_ID_info_RodCones$Type == 'RodCones'),]$Total = sum(COUNT_ALL_ID_info_RodCones[(COUNT_ALL_ID_info_RodCones$Species_ID_ID == ID) & (COUNT_ALL_ID_info_RodCones$Type == 'RodCones'),]$values_plot)
  COUNT_ALL_ID_info_RodCones$Percent = (COUNT_ALL_ID_info_RodCones$values_plot / COUNT_ALL_ID_info_RodCones$Total)*100
  
}

Opsins_tmp_RodCones = unique(subset(COUNT_ALL_ID_info_RodCones, select=c(Opsin,Colors)))
Opsins_final_RodCones = left_join(data.frame(Opsin=c('SWS1', 'SWS2B', 'SWS2A', 'RH2B', 'RH2Ab', 'RH2Aa', 'LWS', 'RH1')), Opsins_tmp_RodCones)
Opsins_order_RodCones = Opsins_final_RodCones$Opsin
Opsins_cols_order_RodCones = Opsins_final_RodCones$Colors


##########################################################################################
# Single and double cone opsins together
##########################################################################################

COUNT_ALL_ID_info_Cones_SingleDouble = COUNT_ALL_ID_info_Cones
COUNT_ALL_ID_info_Cones_SingleDouble[COUNT_ALL_ID_info_Cones_SingleDouble$Opsin %in% c('SWS1', 'SWS2B', 'SWS2A'),]$Type = 'SingleCones'
COUNT_ALL_ID_info_Cones_SingleDouble[COUNT_ALL_ID_info_Cones_SingleDouble$Opsin %in% c('RH2B', 'RH2Ab', 'RH2Aa', 'LWS'),]$Type = 'DoubleCones'

for (ID in unique(Method_count$Species_ID_ID)){
  COUNT_ALL_ID_info_Cones_SingleDouble[(COUNT_ALL_ID_info_Cones_SingleDouble$Species_ID_ID == ID) & (COUNT_ALL_ID_info_Cones_SingleDouble$Type == 'SingleCones'),]$Total = sum(COUNT_ALL_ID_info_Cones_SingleDouble[(COUNT_ALL_ID_info_Cones_SingleDouble$Species_ID_ID == ID) & (COUNT_ALL_ID_info_Cones_SingleDouble$Type == 'SingleCones'),]$values_plot)
  COUNT_ALL_ID_info_Cones_SingleDouble[(COUNT_ALL_ID_info_Cones_SingleDouble$Species_ID_ID == ID) & (COUNT_ALL_ID_info_Cones_SingleDouble$Type == 'DoubleCones'),]$Total = sum(COUNT_ALL_ID_info_Cones_SingleDouble[(COUNT_ALL_ID_info_Cones_SingleDouble$Species_ID_ID == ID) & (COUNT_ALL_ID_info_Cones_SingleDouble$Type == 'DoubleCones'),]$values_plot)
  
  COUNT_ALL_ID_info_Cones_SingleDouble$Percent = (COUNT_ALL_ID_info_Cones_SingleDouble$values_plot / COUNT_ALL_ID_info_Cones_SingleDouble$Total)*100
  
}

COUNT_ALL_ID_info_Cones_SingleDouble[is.na(COUNT_ALL_ID_info_Cones_SingleDouble$Percent),]$Percent = 0


##########################################################################################
# Save variables and RData
##########################################################################################
assign(paste0(Method, '_', 'Method_count', '_exons'), Method_count)
assign(paste0(Method, '_', 'ALL_ID_filter', '_exons'), ALL_ID_filter)
assign(paste0(Method, '_', 'ALL_ID_TPM', '_exons'), ALL_ID_TPM)
assign(paste0(Method, '_', 'ALL_ID_init', '_exons'), ALL_ID_init)
assign(paste0(Method, '_', 'ALL_ID_stats', '_exons'), ALL_ID_stats)
############################################

# Remove variables
rm(Method_count, ALL_ID_filter, COUNT_ALL_ID, COUNT_ALL_ID_info, COUNT_ALL_ID_info_Cones, COUNT_ALL_ID_info_RodCones, COUNT_ALL_ID_info_Cones_SingleDouble,
   ALL_ID_init, ALL_ID_stats,
   ALL_ID_filter_T, ALL_ID_TPM, ALL_ID_TPM_ID)

# Save RData
save.image(file = paste0(PATH, Method, '_exons_FilteredForTPMandVST_withRH2As.RData'))

