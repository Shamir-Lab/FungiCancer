#https://bioconductor.org/packages/release/bioc/vignettes/MMUPHin/inst/doc/MMUPHin.html

#metadata <- read.delim("metadata_fungi_14495samples.tsv", sep = "\t", header = TRUE, fill = TRUE, row.names = 1)
# Load transposed raw counts data
#raw_counts <- read.csv("fungi_rawcounts_only_pt.csv", header = TRUE, row.names = 1)
#mismatched_samples <- setdiff(rownames(metadata), rownames(raw_counts))
#metadata <- metadata[-which(rownames(metadata) %in% mismatched_samples), ]
#t_raw_counts = t(raw_counts)
#fit_adjust_batch <- adjust_batch(feature_abd = t_raw_counts,
#                                 batch = c("data_submitting_center_label"),
#                                 data = metadata,
#                                 control = list(verbose = TRUE))
#
# CRC_abd_adj <- fit_adjust_batch$feature_abd_adj


#write.csv(t(CRC_abd_adj), "mmuphin_adjusted_fungi_2022_onlypt.csv", row.names = TRUE, quote=FALSE)


###############################################
#if (!"BiocManager" %in% rownames(installed.packages()))
#  install.packages("BiocManager")
# BiocManager::install("MBECS")
#BiocManager::install("mixOmics")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("MMUPHin")


library(magrittr)
library(dplyr)
library(mixOmics)
library(ggplot2)
library(MMUPHin) # adjust_batch
library(MBECS) #bmc
library(compositions)
library(Tjazi)  # clr_c
#library(mixOmics) 
library(ggfortify)
library(sva) #combat_seq
library(limma) #removeBatchEffect
library(magrittr)
library(dplyr)
library(ggplot2)
library(vegan)
library(rlang)
library(edgeR)
library(PLSDAbatch)

batch_feature <- 'data_submitting_center_label'
is_permanova_calc <- TRUE
specimen_type <- "RNA-Seq"#"WGS" # #
study <- 'NH22' #'SALZ' #'NH22' #
cohort_defintion <- 'no_harvard' #'all' #'no_harvard' #'iptw' # 'all' 
  
analysis_path <-  "C:\\Users\\danco\\My Drive\\Master\\Projects\\Fungi_Cancer (Thomy Margalit)\\code\\"
setwd(analysis_path)

# Load data
if (study == 'NH22'){
  metadata <- read.delim("metadata_fungi_14495samples.tsv", sep = "\t", header = TRUE, fill = TRUE, row.names = 1)
}
if (study == 'SALZ'){
  metadata <- read.csv("metadata_for_WGS.csv", header = TRUE, row.names = 1)
  #metadata_names <- readxl::read_excel("TableS12_Metadata-TCGA-WGS-5734-Samples.xlsx")
}

# Load transposed raw counts data
if (study == 'NH22'){
  raw_counts <- read.csv("count_data_fungi_decontaminated_raw.csv", header = TRUE, row.names = 1)# fungi_rawcounts_only_pt.csv
}
if (study == 'SALZ'){
  raw_counts <- readxl::read_excel("TableS4_Kraken-TCGA-WGS-5734-Samples-557-Species.Fungi_RefSeq.RawCounts.xlsx")
}

#raw_counts <- read.csv("raw_fungi_for_bc.csv", header = TRUE, row.names = 1)# fungi_rawcounts_only_pt.csv
#mismatched_samples <- setdiff(rownames(metadata), rownames(raw_counts))
#metadata <- metadata[-which(rownames(metadata) %in% mismatched_samples), ]

# filter only black and white data
#metadata_race <- metadata[metadata$race %in% c("BLACK OR AFRICAN AMERICAN","WHITE","ASIAN"),]
#raw_count.race <- raw_counts[rownames(raw_counts) %in% rownames(metadata_race),]

# filter only RNA-Seq
if (study == 'NH22'){
  metadata.race.rna <- metadata[metadata$experimental_strategy %in% c(specimen_type) & metadata$sample_type %in% c("Primary Tumor"),]
  raw_count.race.rna <- raw_counts[rownames(raw_counts) %in% rownames(metadata.race.rna),]
}
if (study == 'SALZ'){
  metadata.race.rna <- metadata[metadata$sample_type %in% c("primary_tumor"),]
  raw_count.race.rna <- raw_counts[raw_counts$HopkinsID %in% metadata.race.rna$HopkinsID,]
  rownames_stored  <- raw_count.race.rna$HopkinsID
  raw_count.race.rna <- raw_count.race.rna %>% dplyr::select(-HopkinsID)
  rownames(raw_count.race.rna) <- rownames_stored
  metadata.race.rna <- metadata.race.rna %>% rename(data_submitting_center_label=data_submitting_center, 
                                                    sample_name=sample_id,
                                                    gender = sex,
                                                    age_at_diagnosis = age)
}
metadata.race.rna$race <- factor(metadata.race.rna$race)

write.csv(metadata.race.rna[,c('sample_name','data_submitting_center_label')], paste('batch_label_',specimen_type,'_',cohort_defintion,'_',study,'.csv',sep=''))

###
# Managing batch effects in microbiome data 
# https://drive.google.com/file/d/1pWLa4UYTfOOyb2IBt49gNHWuttRXcrc6/view
# https://evayiwenwang.github.io/Managing_batch_effects/

##############################
##### Data preprocessing #####
##############################

# (1) Prefilter the count data to remove features with excess zeroes across all samples
# We use a prefiltering step to remove OTUs for which the sum of counts are below a set threshold (0.001%)
# compared to the total sum of all counts (Arumugam et al. 2011).

#Autoplot - https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html

print(dim(raw_count.race.rna))

min_per =  0 # 0.001
race.rna.index.keep <- which(colSums(raw_count.race.rna)*100/(sum(colSums(raw_count.race.rna))) > min_per)
if (study == 'NH22'){
  raw_count.race.rna.keep <- raw_count.race.rna[, race.rna.index.keep]
}
if (study == 'SALZ'){
  raw_count.race.rna[,'temp'] <- rownames(raw_count.race.rna)
  raw_count.race.rna.keep <- raw_count.race.rna[, c(names(race.rna.index.keep),'temp')]
  filtered_rows <- raw_count.race.rna.keep[,'temp']
  raw_count.race.rna.keep <- raw_count.race.rna.keep %>% dplyr::select(-temp)
  rownames(raw_count.race.rna.keep) <- filtered_rows$temp
}


if (cohort_defintion == 'all'){
} else {
  if(cohort_defintion == 'iptw'){
    if (study == 'NH22'){
      #take only subset of data that has the confounders of IPTW
      df_iptw <- read.csv(paste('raw_fungi_for_analysis_',specimen_type,'.csv',sep=''))
      raw_count.race.rna.keep <- raw_count.race.rna.keep[rownames(raw_count.race.rna.keep) %in% df_iptw$Sample,]
      metadata.race.rna <- metadata.race.rna[rownames(metadata.race.rna) %in% df_iptw$Sample,]
    } 
    if (study == 'SALZ'){
      iptw_cols <- c("gender", "age_at_diagnosis", "race","histological_diagnosis_label","pathologic_stage_label")
      df_iptw <- (metadata.race.rna[c(iptw_cols,c("sample_name","HopkinsID","data_submitting_center_label"))])
      df_iptw[df_iptw$pathologic_stage_label == '','pathologic_stage_label'] <- NA
      df_iptw[df_iptw$data_submitting_center_label == '','data_submitting_center_label'] <- NA
      df_iptw[df_iptw$race == 'not reported','race'] <- NA
      df_iptw <- df_iptw[complete.cases(df_iptw[, c(iptw_cols,c(batch_feature))]),]
      
      raw_count.race.rna.keep[,'temp'] <- rownames(raw_count.race.rna.keep)
      raw_count.race.rna.keep <- raw_count.race.rna.keep[raw_count.race.rna.keep$temp %in% df_iptw$HopkinsID,]
      filtered_rows <- raw_count.race.rna.keep[,'temp']
      raw_count.race.rna.keep <- raw_count.race.rna.keep %>% dplyr::select(-temp)
      rownames(raw_count.race.rna.keep) <- filtered_rows$temp
      metadata.race.rna <- metadata.race.rna[metadata.race.rna$HopkinsID %in% df_iptw$HopkinsID,]
      rownames(metadata.race.rna) <- metadata.race.rna$HopkinsID
    }

  } else{
    #exclude samples from Broad Institute of MIT and Harvard
    harvard_ids <- rownames(metadata.race.rna[metadata.race.rna$data_submitting_center_label == 'Broad Institute of MIT and Harvard',])
    raw_count.race.rna.keep <- raw_count.race.rna.keep[!rownames(raw_count.race.rna.keep) %in% harvard_ids,]
    metadata.race.rna <- metadata.race.rna[!rownames(metadata.race.rna) %in% harvard_ids,]
  }
}


#exclude samples with total counts = zero
total_row_sum <- rowSums(raw_count.race.rna.keep)
raw_count.race.rna.keep[,'temp'] <- rownames(raw_count.race.rna.keep)
raw_count.race.rna.keep <- raw_count.race.rna.keep[total_row_sum != 0,]
filtered_rows <- raw_count.race.rna.keep[,'temp']
raw_count.race.rna.keep <- raw_count.race.rna.keep %>% dplyr::select(-temp)
rownames(raw_count.race.rna.keep) <- filtered_rows$temp
metadata.race.rna <- metadata.race.rna[metadata.race.rna$HopkinsID %in% rownames(raw_count.race.rna.keep),]
rownames(metadata.race.rna) <- filtered_rows$temp

#'raw_counts'
res_df <- data.frame(transformation_method = character(), bc_method = character(), permanova_p_val = numeric())

for (method_transform in c('relative_abundance','clr_c','clr_offset','log_cpm_quantile')){
  output_transformed_counts <- transform_counts(raw_count.race.rna.keep, metadata.race.rna, method_transform, is_permanova_calc,batch_feature)
  raw_count.race.rna.keep.transformed <- output_transformed_counts[[1]]
  p_val <- output_transformed_counts[[2]]
  
  res_df <- rbind(res_df, data.frame(transformation_method =  method_transform, bc_method = 'null', permanova_p_val = p_val))
  message(sprintf("transformation method: %s, Permanova p-val: %s", method_transform,p_val))
  
  write.csv(raw_count.race.rna.keep.transformed, paste('adjusted_',specimen_type,'_',cohort_defintion,'_',method_transform,'.csv',sep=''))
  for (method_bc in c('bmc','rbe','combat','mmuphin','plsda')){
    if(!((method_bc == 'plsda') & (method_transform %in% c('relative_abundance','log_cpm_quantile')))){
      message(sprintf("batch correction method: %s", method_bc))
      output_batch_correction <- batch_correction_fun(raw_count.race.rna.keep.transformed,raw_count.race.rna.keep,metadata.race.rna, method_bc, method_transform,batch_feature)
      raw_count.race.rna.keep.transformed.bc <- output_batch_correction[[1]]
      p_val <- output_batch_correction[[2]]
      
      res_df <- rbind(res_df, data.frame(transformation_method =  method_transform, bc_method = method_bc, permanova_p_val = p_val))
      message(sprintf("transformation method: %s, bc method: %s, Permanova p-val: %s", method_transform,method_bc,p_val))
      write.csv(raw_count.race.rna.keep.transformed.bc, paste('adjusted_',specimen_type,'_',cohort_defintion,'_',method_transform,'_',method_bc,'_',study,'.csv',sep=''))
      
    }
  }
  
}

write.csv(res_df, paste('permanova_p_values_no_plsda_',specimen_type,'_',cohort_defintion,'_',method_transform,'_',method_bc,'.csv',sep=''))


method_transform <- 'clr_c'
method_bc <- 'plsda'
p_val_clr_c_plsda <- calc_bc_permanova_p_val(metadata.race.rna,method_transform,method_bc)
res_df <- rbind(res_df, data.frame(transformation_method =  method_transform, bc_method = method_bc, permanova_p_val = p_val_clr_c_plsda))
method_transform <- 'clr_offset'
p_val_clr_offset_plsda <- calc_bc_permanova_p_val(metadata.race.rna,method_transform,method_bc)
res_df <- rbind(res_df, data.frame(transformation_method =  method_transform, bc_method = method_bc, permanova_p_val = p_val_clr_offset_plsda))

write.csv(res_df, 'permanova_p_values.csv')


# PCA
#permanova <- 'null'
#p_val <- 'null'

####################################
###### Apply transformations #######
####################################
# Raw Counts
raw_count.race.rna.keep <- transform_counts(raw_count.race.rna.keep, metadata.race.rna, 'raw_counts', TRUE)
# Relative abundance transformation
raw_count.race.rna.keep.ra <- transform_counts(raw_count.race.rna.keep, metadata.race.rna, 'relative_abundance', TRUE)
# Log-ratio transformation with Centered Log Ratio (CLR). clr_c impute zeroes and perform a centered log-ratio.
#  https://rdrr.io/github/thomazbastiaanssen/Tjazi/src/R/clr_lite.R
raw_count.race.rna.keep.clr_c <- transform_counts(raw_count.race.rna.keep, metadata.race.rna, 'clr_c', TRUE)
# clr with offset (Add an offset of 1 to the whole data matrix — note that this is not ideal but provides a pracitical way to handle zero counts.)
raw_count.race.rna.keep.clr_offset <- transform_counts(raw_count.race.rna.keep, metadata.race.rna, 'clr_offset', TRUE)


transform_counts <- function(counts, metadata.race.rna, method_transform, is_permanova_calc, batch_feature){
  if (method_transform == 'raw_counts'){
    counts <- counts
    permanova <- adonis2(counts ~ eval(rlang::sym(batch_feature)),
                         data = metadata.race.rna, permutations=99, method = "euclidean")
  }
  if (method_transform == 'relative_abundance'){
    counts <- t(apply(counts,1, function(x) x/sum(x)))
    permanova <- adonis2(counts ~ eval(rlang::sym(batch_feature)),
                         data = metadata.race.rna, permutations=99, method = "euclidean")
  }
  if (method_transform == 'clr_c'){
    #https://github.com/thomazbastiaanssen/Tjazi
    counts <- t(clr_c(t(counts)))
    permanova <- adonis2(counts ~ eval(rlang::sym(batch_feature)),
                         data = metadata.race.rna, permutations=99, method = "euclidean")
  }
  if (method_transform == 'clr_offset'){
    offset <- 1
    counts <- counts + offset
    counts <- logratio.transfo(counts, logratio = 'CLR')
    class(counts) <- 'matrix' 
    permanova <- adonis2(counts ~ eval(rlang::sym(batch_feature)),
                         data = metadata.race.rna, permutations=99, method = "euclidean")
  }
  if (method_transform == 'log_cpm_quantile'){
    dge <- DGEList(t(counts)) # DGEList object from a table of counts (rows=features, columns=samples)
    counts <- t(voom(dge,  normalize.method="quantile")$E)
    permanova <- adonis2(counts ~ eval(rlang::sym(batch_feature)),
                         data = metadata.race.rna, permutations=99, method = "euclidean")
  }
  
  if (is_permanova_calc){
    #p_val <- as.data.frame(permanova$aov.tab)["data_submitting_center_label", "Pr(>F)"] #adonis and not adonis2
    p_val <- permanova$`Pr(>F)`[1]
    pca_res <- prcomp(counts, scale. = TRUE)
    pca_plot <-autoplot(pca_res, data =metadata.race.rna, colour = batch_feature, frame = TRUE, frame.type = 'norm',main = paste(method_transform, "\n(Permanova P=", p_val,')',sep=''))+ 
      theme(legend.position="top",legend.text=element_text(size=6),legend.title=element_text(size=8)) + guides(colour = guide_legend(nrow = 3))
    ggsave(paste('pca_plot_',method_transform,'.png',sep=''), pca_plot)
  }
  else{
    p_val <- 99999
  }
  
  return(list(counts, p_val))
}

batch_correction_fun <- function(transformed_counts,counts,metadata.race.rna, method_bc, method_transform, batch_feature){
  if (method_bc == 'bmc'){
    transformed_counts <-  batch_correction_bmc_mbecs(transformed_counts, metadata.race.rna, batch_feature)
    #transformed_counts <- batch_correction_bmc(transformed_counts, metadata.race.rna)
    permanova <- adonis2(transformed_counts ~ eval(rlang::sym(batch_feature)),
                         data = metadata.race.rna, permutations=99, method = "euclidean")
  }
  if (method_bc == 'rbe'){
    transformed_counts <- t(limma::removeBatchEffect(t(transformed_counts), batch = metadata.race.rna$data_submitting_center_label))
    permanova <- adonis2(transformed_counts ~ eval(rlang::sym(batch_feature)),
                         data = metadata.race.rna, permutations=99, method = "euclidean")
  }
  if (method_bc == 'combat'){
    transformed_counts <- t(ComBat(t(transformed_counts), batch=metadata.race.rna$data_submitting_center_label,mod=NULL))
    permanova <- adonis2(transformed_counts ~ eval(rlang::sym(batch_feature)),
                         data = metadata.race.rna, permutations=99, method = "euclidean")
  }
  if (method_bc == 'plsda'){
    if (method_transform %in% c('clr_offset','clr_c')){
      transformed_counts <- batch_correct_plsda(transformed_counts, metadata.race.rna[,batch_feature])
      permanova <- adonis2(transformed_counts ~ eval(rlang::sym(batch_feature)),
                           data = metadata.race.rna, permutations=99, method = "euclidean")
    }
  }
  if (method_bc == 'mmuphin'){
    #https://bioconductor.org/packages/release/bioc/vignettes/MMUPHin/inst/doc/MMUPHin.html
    if (method_transform %in% c('raw_counts','relative_abundance')){
      if(method_transform == 'relative_abundance'){
        counts <- t(apply(counts,1, function(x) x/sum(x)))
      }
      fit_adjust_batch <- adjust_batch(feature_abd = t(counts),
                                       batch = c(batch_feature),
                                       data = metadata.race.rna)
      transformed_counts <- t(fit_adjust_batch$feature_abd_adj)
    }
    if (method_transform %in% c('clr_offset','clr_c')){
      fit_adjust_batch <- adjust_batch(feature_abd = t(counts),
                                       batch = c(batch_feature),
                                       data = metadata.race.rna)
      # transform adjusted data
      transformed_counts <- transform_counts(t(fit_adjust_batch$feature_abd_adj), metadata.race.rna, method_transform, FALSE)
      transformed_counts <- transformed_counts[[1]]
    }
    permanova <- adonis2(transformed_counts ~ eval(rlang::sym(batch_feature)),
                         data = metadata.race.rna, permutations=99, method = "euclidean")
  }
  
  p_val <- permanova$`Pr(>F)`[1]
  pca_res <- prcomp(transformed_counts, scale. = TRUE)
  pca_plot <-autoplot(pca_res, data =metadata.race.rna, colour = batch_feature, frame = TRUE, frame.type = 'norm',main = paste(method_transform,'+',method_bc, "\n(Permanova P=", p_val,')',sep=''))+ 
    theme(legend.position="top",legend.text=element_text(size=6),legend.title=element_text(size=8)) + guides(colour = guide_legend(nrow = 3))
  ggsave(paste('pca_plot_',method_transform,'_',method_bc,'.png',sep=''), pca_plot)
  return(list(transformed_counts,p_val))
}

calc_bc_permanova_p_val <- function(metadata.race.rna,method_transform,method_bc){
  temp <- metadata.race.rna
  temp['X'] <- rownames(temp)
  transformed_counts <- read.csv(paste('adjuste_fungi_onlypt_',method_transform,'_',method_bc,'.csv',sep=''))
  transformed_counts <- merge(transformed_counts, temp[,c('X',"data_submitting_center_label")],  by='X')
  transformed_counts$data_submitting_center_label <- as.numeric(factor(transformed_counts$data_submitting_center_label))
  row.names(transformed_counts) <- transformed_counts[,'X']
  transformed_counts <- transformed_counts[ , -which(names(transformed_counts) %in% c("X"))]
  permanova <- adonis2(transformed_counts ~ eval(rlang::sym(batch_feature)),
                       data = metadata.race.rna, permutations=99, method = "euclidean")
  p_val <- permanova$`Pr(>F)`[1]
  
  pca_res <- prcomp(transformed_counts, scale. = TRUE)
  pca_plot <-autoplot(pca_res, data =metadata.race.rna, colour = batch_feature, frame = TRUE, frame.type = 'norm',main = paste(method_transform,'+',method_bc, "\n(Permanova P=", p_val,')',sep=''))+ 
    theme(legend.position="top",legend.text=element_text(size=6),legend.title=element_text(size=8)) + guides(colour = guide_legend(nrow = 3))
  ggsave(paste('pca_plot_',method_transform,'_',method_bc,'.png',sep=''), pca_plot)
  return (p_val) 
}


#Batch mean centering
batch_correction_bmc <- function(raw_count.race.rna.keep, metadata.race.rna){
  raw_count.race.rna.keep.b1 <- scale(raw_count.race.rna.keep[metadata.race.rna$data_submitting_center_label == 'University of North Carolina', ], center = TRUE, scale = FALSE)
  raw_count.race.rna.keep.b2 <- scale(raw_count.race.rna.keep[metadata.race.rna$data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre", ], center = TRUE, scale = FALSE)
  raw_count.race.rna.keep.b3 <- scale(raw_count.race.rna.keep[metadata.race.rna$data_submitting_center_label == "Broad Institute of MIT and Harvard", ], center = TRUE, scale = FALSE)
  raw_count.race.rna.keep.bmc <- rbind(raw_count.race.rna.keep.b1,raw_count.race.rna.keep.b2,raw_count.race.rna.keep.b3)
  #raw_count.race.rna.keep.clr.bmc <- raw_count.race.rna.keep.clr.bmc[rownames(raw_count.race.rna.keep.clr), ]
  bmc_adjusted <- raw_count.race.rna.keep.bmc[rownames(raw_count.race.rna.keep), ]
  return(bmc_adjusted)
}
#https://yury-zablotski.netlify.app/post/mixed-models/#simple-random-slope-model

# batch mean centering - mbecs
# MBECS  - https://github.com/rmolbrich/MBECS

batch_correction_bmc_mbecs <- function(raw_count.race.rna.keep, metadata.race.rna, batch_feature){
  data("dummy.list")
  new_obj <- dummy.list
  new_obj$cnts <- raw_count.race.rna.keep
  new_obj$meta <- metadata.race.rna
  mbec.obj <- mbecProcessInput(new_obj)
  mbec.obj <- mbecCorrection(mbec.obj, model.vars=c(batch_feature), 
                             method = "bmc", type='otu')
  bmc_adjusted <- t(mbec.obj@corrections$bmc)
  return(bmc_adjusted)
}

# Define the batch correction function for PLSDA
batch_correct_plsda <- function(df, batches, balance=FALSE) {
  result <- PLSDA_batch(df, NULL, batches, balance)
  return(result$X.nobatch)
}

######################################
###### Apply Batch-corrections #######
######################################

raw_count.race.rna.keep.clr_c.bmc  <- batch_correction_fun(raw_count.race.rna.keep.clr_c,raw_count.race.rna.keep,metadata.race.rna, 'bmc', 'clr_c')
# RemoveBatchEffects - #https://rdrr.io/bioc/limma/man/removeBatchEffect.html
raw_count.race.rna.keep.clr_c.rbe  <- batch_correction_fun(raw_count.race.rna.keep.clr_c,raw_count.race.rna.keep,metadata.race.rna, 'rbe', 'clr_c')

#combat
raw_count.race.rna.keep.clr_c.combat  <- batch_correction_fun(raw_count.race.rna.keep.clr_c,raw_count.race.rna.keep,metadata.race.rna, 'combat', 'clr_c')

#mmuphin
# https://rdrr.io/bioc/MMUPHin/man/adjust_batch.html
# https://bioconductor.org/packages/release/bioc/vignettes/MMUPHin/inst/doc/MMUPHin.html#4_Performing_batch_(study)_effect_adjustment_with_adjust_batch
raw_count.race.rna.keep.clr_c.mmuphin  <- batch_correction_fun(raw_count.race.rna.keep.clr_c,raw_count.race.rna.keep,metadata.race.rna, 'mmuphin', 'clr_c')




##################
### COMBAT-seq ###
##################
# https://rdrr.io/bioc/sva/man/ComBat_seq.html
# https://github.com/zhangyuqing/ComBat-seq
#https://academic.oup.com/nargab/article/2/3/lqaa078/5909519
# Add offset of +1 to counts
#offset_param <- 0.5
#combat_adjusted <- sva::ComBat_seq(t(raw_count.race.rna.keep)+offset_param, batch=metadata.race.rna$data_submitting_center_label,group=NULL, full_mod=FALSE)

##################
## alternative BMC
##################

new_obj <- dummy.list
new_obj$cnts <- raw_count.race.rna.keep
new_obj$meta <- metadata.race.rna
mbec.obj <- mbecProcessInput(new_obj)
mbec.obj <- mbecCorrection(mbec.obj, model.vars=c("data_submitting_center_label","race"), 
                           method = "bmc", type='otu')
bmc_adjusted <- t(mbec.obj@corrections$bmc)






metadata_race <- metadata[metadata$race %in% c("BLACK OR AFRICAN AMERICAN","WHITE"),]
raw_count_race <- raw_counts[rownames(raw_counts) %in% rownames(metadata_race),]
race.index.keep <- which(colSums(raw_count_race)*100/(sum(colSums(raw_count_race))) > 0.001)
raw_count_race.keep <- raw_count_race[, race.index.keep]
print(dim(raw_count_race.keep))

# PCA
pca_res <- prcomp(raw_count_race.keep, scale. = TRUE)
autoplot(pca_res, data =metadata_race, colour = 'data_submitting_center_label', frame = TRUE, frame.type = 'norm',main = "Raw Counts")

# (2) Add an offset of 1 to the whole data matrix — note that this is not ideal but provides a pracitical way to handle zero counts.

## TBD

# (3) Log-ratio transformation with Centered Log Ratio (CLR). clr_c impute zeroes and perform a centered log-ratio.
transformed_raw_count_race <- clr_c(raw_count_race.keep)

rbe_adjusted_two_batches <- removeBatchEffect(t(transformed_raw_count_race), batch = metadata_race$data_submitting_center_label)

pca_res <- prcomp(t(rbe_adjusted), scale. = TRUE)
autoplot(pca_res, data =metadata_race, colour = 'data_submitting_center_label', frame = TRUE, frame.type = 'norm',main = "removeBatchEffect")



mbecReportPrelim(input.obj=mbec.obj, model.vars=c("data_submitting_center_label","race"), 
                 type="clr")

#mbec.obj <- mbecCorrection(mbec.obj, model.vars=c("data_submitting_center_label","race"), 
#                                     method = "bat", type = "clr")

mbec.obj <- mbecCorrection(mbec.obj, model.vars=c("data_submitting_center_label","race"), 
                           method = "rbe", type='otu')

#mbec.obj <- mbecRunCorrections(mbec.obj, model.vars=c("data_submitting_center_label","race"),
#                               method=c("ruv3","rbe","bmc","pn","svd"), 
#                               type = "clr")
new_obj <- dummy.list
new_obj$cnts <- raw_count_race
metadata_race$race <- factor(metadata_race$race)
new_obj$meta <- metadata_race
mbec.obj_2 <- mbecProcessInput(new_obj)

mbec.obj_2 <- mbecCorrection(mbec.obj_2, model.vars=c("data_submitting_center_label","race"), 
                           method = "pn", type='otu')
#Apply PCA
library(ggfortify)

pca_res <- prcomp(raw_count_race, scale. = TRUE)
autoplot(pca_res, data =metadata_race, colour = 'data_submitting_center_label')


pca_res <- prcomp(transformed_raw_count_race, scale. = TRUE)
autoplot(pca_res, data =metadata_race, colour = 'data_submitting_center_label')


pca_res <- prcomp(t(mbec.obj@corrections$rbe), scale. = TRUE)
autoplot(pca_res, data =metadata_race, colour = 'data_submitting_center_label')


pca_res <- prcomp(t(mmuphin_adjusted), scale. = TRUE)
autoplot(pca_res, data =metadata_race, colour = 'data_submitting_center_label')

pca_res <- prcomp(t(mbec.obj_2@corrections$pn), scale. = TRUE)
autoplot(pca_res, data =metadata_race, colour = 'data_submitting_center_label')

#####################
### PLDSA - Batch ###
#####################
# https://rdrr.io/github/EvaYiwenWang/PLSDAbatch/
# https://support.bioconductor.org/p/133791/#133825
# https://www.biorxiv.org/content/10.1101/2020.10.27.358283v1.full



