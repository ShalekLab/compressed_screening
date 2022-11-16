
#### title: "QC filtering for PDAC compressed screen plate1"
#### author: "BEM"
  ### Metadata
  
### Metadata
### From EA3.16 PDAC orgnaoid compressed ligand screen
### each sample contains 12 wells - hashed via BD antibodies

# 
# ### Preprocessing
# 
# Run Basic QC (per array) on the raw Cell X gene matricies

Sys.setenv(PATH = paste("/Users/bemead/anaconda3/bin", Sys.getenv("PATH"), sep=":"))
setwd("~/Cloud/Dropbox (MIT)/Shalek Lab Projects/MITGSK/Aim3-Screening/EXPT_A3.16_BEM_10.2021/3_Analysis/")

source("1_Preprocess/Preprocess_QC_BEM_v1.5_HASHING.R")

sample_sheet = "1_Preprocess/SSwl_q1_array_list.txt"
figdir = "1_Preprocess/plots/starsolo_wl_g500"

Objs = QC_preprocess_filteronly(data_type = 'StarSolo',
                                max_mito = 50, min_feat = 500, 
                                max_UMI = 30000, min_UMI =500,
                                sample_sheet = sample_sheet, figdir = figdir,
                                drop_unanno = T)

# Export 1st round of objects - 1st plate (q1)
saveRDS(Objs, '1_Preprocess/PDAC_CS_SSwl_q1_QC_prepoc_g500_BEM_20220901.rds', compress= T)
