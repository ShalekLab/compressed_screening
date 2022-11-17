#### Created By Ben Mead -- November 2019
### Last update 05.31.2022
### update for processing cell hashing seq-well expts
### V1.5

#### QC PIPELINE FUNCTIONS FOR Honeycomb MATRICIES ALIGNED ON TERRA

Sys.setenv(PATH = paste("/Users/bemead/anaconda3/bin", Sys.getenv("PATH"), sep=":"))
require(tidyverse, quietly = T)
require(Seurat, quietly = T)
require(DoubletFinder, quietly = T)
require(fields, quietly = T)
require(KernSmooth, quietly = T)
require(modes, quietly = T)
require(ROCR, quietly = T)
require(rhdf5, quietly = T)
require(SparseM, quietly = T)
require(deMULTIplex, quietly = T)

# quiet output
hush = function(code){
  sink("/dev/null")
  tmp = code
  sink()
  return(tmp)
}

#### QC_preprocess fn (QC filter only)
QC_preprocess_filteronly = function(sample_sheet = "1_Preprocess/array_list.txt", data_type = 'DropSeqTools', 
                         min_feat = 500, min_UMI = 50, max_UMI = 30000, max_mito = 50, figdir = "1_Preprocess/plots", drop_unanno = F){
  
  # Read in array file
  samples = read.table(sample_sheet, header = T)
  n_sample = length(samples$array)
  
  # Run initial QC filtering
  names = c()
  Objs = list()
  for (i in 1:n_sample){
    array = as.character(samples$array[i])
    dir = as.character(samples$dir[i])
    
    # Build  object
    if (data_type == 'DropSeqTools'){
      data = read.table(dir, header = T, row.names = 1)
    } else if (data_type == 'StarSolo'){
      data = Read10X(dir)
    } else if (data_type == 'CellBender'){
      data.shape = h5read(file = dir, name = "background_removed/shape")
      data.genes = h5read(file = dir, name = "background_removed/gene_names") %>% as.character()
      data.barcodes = h5read(file = dir, name = "background_removed/barcodes") %>% as.character()
      data.data = h5read(file = dir, name = "background_removed/data")
      data.indices = h5read(file = dir, name = "background_removed/indices")
      data.indptr = h5read(file = dir, name = "background_removed/indptr")
      
      data = new("matrix.csc", ra = as.numeric(data.data), ja = as.integer(1+data.indices), ia = as.integer(1+data.indptr), dimension = as.integer(data.shape)) %>% as("dgCMatrix")
      rownames(data) = data.genes
      colnames(data) = data.barcodes
    }
    
    # Filter unannotated
    if (drop_unanno == T){
      unannotated = grep(pattern = "^A[A-Z][0-9]*.[0-9]", x = rownames(data), value = TRUE)
      LINC = grep(pattern = "^LINC[0-9]*", x = rownames(data), value = TRUE)
      
      data = data[!(rownames(data) %in% unannotated), ]
      data = data[!(rownames(data) %in% LINC), ]
    }
    
    Obj = CreateSeuratObject(data, project = array, min.cells = 5)
    
    # Plot QC
    Obj$percent.mt = PercentageFeatureSet(Obj, pattern = "^MT-")
    Obj$percent.ribo = PercentageFeatureSet(Obj, pattern = "^RP")
    
    
    p4 = FeatureScatter(Obj, 'nCount_RNA', 'nFeature_RNA') +
      geom_vline(aes(xintercept=min_UMI), color='black', linetype=2) +
      geom_vline(aes(xintercept=max_UMI), color='black', linetype=2) +
      geom_hline(aes(yintercept=min_feat), color='black', linetype=2) +
      scale_x_continuous(trans='log10')
    p5 = FeatureScatter(Obj, 'nCount_RNA', 'percent.mt') +
      geom_vline(aes(xintercept=min_UMI), color='black', linetype=2) +
      geom_vline(aes(xintercept=max_UMI), color='black', linetype=2) +
      geom_hline(aes(yintercept=max_mito), color='black', linetype=2) +
      scale_x_continuous(trans='log10')
    p6 = FeatureScatter(Obj, 'nCount_RNA', 'percent.ribo') +
      geom_vline(aes(xintercept=min_UMI), color='black', linetype=2) +
      geom_vline(aes(xintercept=max_UMI), color='black', linetype=2) +
      scale_x_continuous(trans='log10')
    p7 = FeatureScatter(Obj, 'nFeature_RNA', 'percent.mt') +
      geom_vline(aes(xintercept=min_feat), color='black', linetype=2) +
      geom_hline(aes(yintercept=max_mito), color='black', linetype=2)
    p8 = FeatureScatter(Obj, 'nFeature_RNA', 'percent.ribo') +
      geom_vline(aes(xintercept=min_feat), color='black', linetype=2)
    
    # QC filter
    cells.use = Cells(Obj)[which(Obj$nFeature_RNA > min_feat)]
    Obj = subset(Obj, cells = cells.use)
    cells.use = Cells(Obj)[which(Obj$nCount_RNA < max_UMI)]
    Obj = subset(Obj, cells = cells.use)
    cells.use = Cells(Obj)[which(Obj$nCount_RNA > min_UMI)]
    Obj = subset(Obj, cells = cells.use)
    cells.use = Cells(Obj)[which(Obj$percent.mt < max_mito)]
    Obj = subset(Obj, cells = cells.use)
    
    p3 = VlnPlot(Obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
                 ncol = 2, pt.size = 0.5)
    
    # Plot results 
    pdf(paste(figdir,"/",array,"_1.QC.pdf", sep = ""))
    plot(p4)
    plot(p5)
    plot(p6)
    plot(p7)
    plot(p8)
    plot(p3)
    dev.off()
    
    ### MODIFY HERE BASED ON DESIRED METADATA
    # Add relevant metadata
    Obj$array = array
    Obj$percent.mt = PercentageFeatureSet(Obj, pattern = "^MT-")
    Obj$percent.ribo = PercentageFeatureSet(Obj, pattern = "^RP")
    
    names[i] = array
    Objs[i] = Obj
    
  }
  names(Objs) = names
  
  return(Objs)
}

######## DeMux functions

HTODemux_q_scan = function(Obj, array, figdir, n_iter = 50, run_90 = F){
  
  # iterate over  quantiles to find max singlet
  scan_table = data.frame()
  n <- 0
  for (q in (1-exp(seq(log(1e-6), log(0.15), length.out = n_iter)))) {
    
    Obj = HTODemux(Obj, assay = "HTO",
                   positive.quantile = q,
                   kfunc = "kmeans",
                   verbose = F)
    
    class = table(Obj$HTO_classification.global)
    class$q = q
    
    scan_table = bind_rows(scan_table, class)
  }
  
  if (run_90){
    singlet_90th = quantile(scan_table$Singlet, 0.9, na.rm = T)
    max_q = max(scan_table[scan_table$Singlet>=singlet_90th,'q'])
  } else{
    max_q = tail(scan_table[order(scan_table$Singlet), ], n=1)$q  
  }
  
  scan_table = gather(scan_table, 'Doublet', 'Negative', 'Singlet', key = 'class', value = 'cells')
  
  # plot
  p = ggplot(data=scan_table, aes(x=log(1-q), y=cells, color=class)) + geom_line()+
    geom_vline(xintercept=log(1-max_q), lty=2) + scale_color_manual(values=c("red","black","blue"))
  
  pdf(paste(figdir,"/",array,"_2.hash_scan.pdf", sep = ""))
  plot(p)
  dev.off()

  # run demux with max q
  Obj = HTODemux(Obj, assay = "HTO",
                 positive.quantile = max_q,
                 kfunc = "kmeans",
                 verbose = F)
  
  return(Obj)
}

#### make some nice hash plots
hash_plots = function(Obj, n_bc=12, array, figdir){
  
   DefaultAssay(Obj) = "HTO"
  
   p1 = RidgePlot(Obj, features = rownames(Obj[["HTO"]])[1:n_bc], 
            group.by = 'hash.ID', stack =T) + NoLegend()
  

  Obj = ScaleData(Obj, features = rownames(Obj), verbose = FALSE)
  Obj = RunPCA(Obj, features = rownames(Obj), approx = FALSE)
  Obj = RunUMAP(Obj, dims = 1:n_bc)
  p2 = DimPlot(Obj, label = F, group.by = 'hash.ID')
  
  p3 = VlnPlot(Obj, features = c("nCount_HTO", "nCount_RNA"), group.by = 'hash.ID' )
  p4 = VlnPlot(Obj, features = c("nCount_HTO", "nCount_RNA"), group.by = 'HTO_classification.global' )
  
  pdf(paste(figdir,"/",array,"_4.hash_class.pdf", sep = ""))
  plot(p1)
  plot(p2)
  plot(p3)
  plot(p4)
  dev.off()
  
}

