
#### title: "Cell hasing demultiplexing CS PDAC plate 1"
#### author: "BEM"
### Metadata

### Metadata
### From EA3.16 PDAC orgnaoid pilot
### From EA3.16 PDAC orgnaoid compressed ligand screen
### each sample contains 12 wells - hashed via BD antibodies

# 
# ### Preprocessing
# 

Sys.setenv(PATH = paste("/Users/bemead/anaconda3/bin", Sys.getenv("PATH"), sep=":"))
setwd("~/Cloud/Dropbox (MIT)/Shalek Lab Projects/MITGSK/Aim3-Screening/EXPT_A3.16_BEM_10.2021/3_Analysis/")

source("1_Preprocess/Preprocess_QC_BEM_v1.5_HASHING.R")

sample_sheet = "1_Preprocess/SSwl_q1_array_list.txt"
figdir = "1_Preprocess/plots/starsolo_wl_g500/"

Objs = readRDS('1_Preprocess/PDAC_CS_SSwl_q1_QC_prepoc_g500_BEM_20220901.rds')

# Load in Hashing libraries
htos = read.table("1_Preprocess/SSwl_q1_hto_list.txt", header = T)
n_sample = length(Objs)

# Run load HTOs into arrays
for (i in 1:n_sample){
  
  dir = as.character(htos$dir[i])
  array = as.character(htos$array[i])
  sample.hto = Read10X(dir, gene.column=1)
  
  # drop unmapped read sample & rename
  sample.hto = sample.hto[-dim(sample.hto)[1],]
  rownames(sample.hto) = gsub("-.*","",rownames(sample.hto))
  
  # drop BCs where HTO nUMI < min_UMI (based on CR 2.2 cell filter)
  colSort = function(data, ...) sapply(data, sort, ...)
  min_UMI = round(quantile(colSort(colSums(sample.hto), decreasing = T)[1:3000], 0.9)/10)
  hto.pass = colnames(sample.hto[,colSums(sample.hto)>min_UMI])
  
  # find matched barcodes
  joint.bcs = intersect(Cells(Objs[[i]]), hto.pass)
  
  # Subset RNA and HTO counts by joint cell barcodes
  Objs[[i]] = subset(Objs[[i]], cells = joint.bcs)
  sample.hto = as.matrix(sample.hto[, joint.bcs])
  
  # Add to seurat object
  Objs[[i]][["HTO"]] = CreateAssayObject(counts = sample.hto)
  
  # normalize HTO
  Objs[[i]] = NormalizeData(Objs[[i]], assay = "HTO", normalization.method = "CLR")
  # drop barcodes where SNR < 1.25
  max_sort = colSort(as.data.frame(GetAssayData(Objs[[i]], assay = "HTO")), decreasing = T)
  SNR = max_sort[1,]/max_sort[2,]
  Objs[[i]] = subset(Objs[[i]], cells = names(SNR)[SNR>1.25])
  # re-normalize
  Objs[[i]] = NormalizeData(Objs[[i]], assay = "HTO", normalization.method = "CLR")
  
  # classify barcodes
  Objs[[i]] = HTODemux_q_scan(Objs[[i]], array, figdir, n_iter =25, run_90=F)
 
  # Specify row
  rows = c('A','B','C','D','E','F','G','H')
  Objs[[i]]$plate_row = rows[i]
   
  }
  
# Merge each array into final object
Obj = merge(Objs[[1]], Objs[2:n_sample])

# make summary plots
levels = c("col-1", "col-2", "col-3", "col-4", "col-5", "col-6", 
           "col-7", "col-8", "col-9", "col-10", "col-11", "col-12", 
           "Doublet", "Negative")
Obj$hash.ID = factor(Obj$hash.ID, levels = levels)
hash_plots(Obj, n_bc=12, 'Merged_Q1', figdir)

#retain only singlets
Obj$plate_column = gsub("...-","",Obj$hash.ID)

Obj = subset(Obj, subset = HTO_classification.global == 'Singlet')
Obj$sample_well = paste(Obj$plate_row,Obj$plate_column, sep="")

# Drop RNA genes where nCell < 10
merge.RNA = GetAssayData(Obj, assay = "RNA")
merge.RNA = merge.RNA[rowSums(merge.RNA)>10,]
Obj[["RNA"]] = CreateAssayObject(counts = merge.RNA)
Obj$percent.mt = PercentageFeatureSet(Obj, pattern = "^MT-")
Obj$percent.ribo = PercentageFeatureSet(Obj, pattern = "^RP")

saveRDS(Obj, '1_Preprocess/PDAC_CS_SSwl_q1_demulti_g500_BEM_20220901.rds', compress= T)
