# library the prerequisite packages
library(Seurat)
library(hdf5r)
library(Matrix)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Write the Hdf5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Instruction
# -- data - An object whose class is Seurat
# -- file - A file path
# -- unmodified.dataset - User defined slots that are not changed to 
# -- ------------------ - improve write-in efficiency
# Function 
# -- Creates a new hdf5 file or Opens an existing one for read/write. 

# however the hdf5 do not support for reclaiming the space
# write the hdf5
# matrix to h5
# matrix to h5
matrix_to_h5 <- function(mat, h5_gp, gp_name = NULL){
  if(!gp_name %in% names(h5_gp)){
    h5mat = h5_gp$create_group(gp_name)
  }
  else{
    h5mat = h5_gp[[gp_name]]
  }
  if('dgCMatrix' %in% class(mat)){
    h5mat[['values']] <- mat@x
    h5mat[['indices']] <- mat@i
    h5mat[['indptr']] <- mat@p
    h5mat[['dims']] <- rev(mat@Dim)
    h5attr(h5mat, 'datatype') <- 'SparseMatrix'
  }
  else if('matrix' %in% class(mat)){
    h5mat[['matrix']] <- mat
    h5mat[['dims']] <- rev(dim(mat))
    h5attr(h5mat, 'datatype') <- 'Array'
  }
}


# dataframe to h5
df_to_h5 <- function(df, h5_anno, anno_dataset = NULL, anno_gp_name = NULL, anno_gp_dataset = NULL){
  # remove the factor
  cate_list <- list()
  for(k in names(df)){
    if(is.factor(df[[k]])){
      obs_cate <- levels(df[[k]])
      df[[k]] <- as.numeric(df[[k]]) - 1 # for 0 begin
      cate_list[[k]] <- obs_cate
    }
  }
  # to h5
  if(is.null(anno_gp_name)){
    h5_anno[[anno_dataset]] <- df
    if(length(cate_list)>0){
      for(j in names(cate_list)){
        h5attr(h5_anno[[anno_dataset]], j) <- cate_list[[j]]
      }
    }
  }
  else{
    if(!anno_gp_name %in% names(h5_anno)){
      anno_gp <- h5_anno$create_group(anno_gp_name)
    }
    else{
      anno_gp <- h5_anno[[anno_gp_name]]
    }
    anno_gp[[anno_gp_dataset]] <- df
    if(length(cate_list)>0){
      for(m in names(cate_list)){
        h5attr(anno_gp[[anno_gp_dataset]], m) <- cate_list[[m]]
      }
    }
  }
}


adata_to_h5 <- function(adata=NULL, h5=NULL){
  # rdata
  rdata <- Seurat::GetAssayData(object = adata, slot = 'counts', assay = 'RNA')
  matrix_to_h5(mat = rdata, h5_gp = h5, gp_name = 'rawData')
  # ndata
  ndata <- Seurat::GetAssayData(object = adata, slot = 'data', assay = 'RNA')
  matrix_to_h5(mat = ndata, h5_gp = h5, gp_name = 'normData')
  # sdata
  sdata <- Seurat::GetAssayData(object = adata, slot = 'scale.data', assay = 'RNA')
  matrix_to_h5(mat = sdata, h5_gp = h5, gp_name = 'scaleData')
  # anno
  annotation <- h5$create_group('annotation')
  # -- observes
  obs = adata@meta.data
  obs[['index']] <- rownames(obs)
  df_to_h5(df = obs, h5_anno = annotation, anno_dataset = 'observes')
  # -- variables
  ndvar = adata@assays$RNA@meta.features
  ndvar[['index']] <- rownames(ndvar)
  hvg <- adata@assays$RNA@var.features
  ndvar[['highly_variable']] <- ndvar[['index']] %in% hvg
  sdvar <- ndvar[hvg, ]
  df_to_h5(df = ndvar, h5_anno = annotation, anno_dataset = 'ndVariables')
  df_to_h5(df = sdvar, h5_anno = annotation, anno_dataset = 'sdVariables')
  if(length(adata@reductions)>0){
    dimReduction <- h5$create_group('dimReduction')
    for(d in names(adata@reductions)){
      D = toupper(d)
      dimReduction[[D]] <- t(adata@reductions[[d]]@cell.embeddings)
    }
  }
  if(length(adata@graphs)>0){
    gra_list <- list(RNA_nn = 'knn', RNA_snn = 'snn')
    graphs <- h5$create_group('graphs')
    graph_df <- adata@graphs
    for(g in names(graph_df)){
      matrix_to_h5(mat=graph_df[[g]], h5_gp = graphs, gp_name = gra_list[[g]])
    }
  }
  if(length(adata@misc)>0){
    meta <- h5$create_group('metadata')
    colo <- grep('color', names(adata@misc), value = TRUE)
    if(length(colo)>0){
      color <- meta$create_group('colors')
      for(co in colo){
        color[[co]] <- adata@misc[[co]]
      }
    }
  }
}

# R write hdf5
R_write_hdf5 <- function(adata = NULL, file = NULL){
  if(is.null(file)){
    stop('No such file or directory') 
  }
  if(class(adata) != 'Seurat'){
    stop('oject ', substitute(adata), ' class is not Seurat object')
  }
  h5 <- H5File$new(filename = file, mode = 'w')
  tryCatch({
    adata_to_h5(adata = adata, h5 = h5)
  }, error = function(err) {
    cat("error!", err, "\n")
  }, finally = {
    h5$close_all()
  }
  )
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read the Hdf5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# 
h5_to_matrix <-  function(gp_name, obs_names, var_names){
  if(!all(rev(gp_name[['dims']][]) == c(length(var_names), length(obs_names)))){
    return(warning('Matrix(SparseMatrix) does not correspond to the annotation dimension, do not read the Matrix', '\n'))
  }else{
    if(h5attr(gp_name, 'datatype') == 'SparseMatrix'){
      mat <- Matrix::sparseMatrix(i = gp_name[['indices']][], 
                                  p = gp_name[['indptr']][], 
                                  x = gp_name[['values']][],
                                  dims = rev(gp_name[['dims']][]), 
                                  index1 = FALSE, 
                                  dimnames = list(var_names, obs_names))
    }
    else if(h5attr(gp_name, 'datatype') == 'Array'){
      if(all(gp_name[['matrix']]$dims == rev(gp_name[['dims']][]))){
        mat <- gp_name[['matrix']][,]
      }
      else{
        mat <- gp_name[['matrix']][,]
        mat <- t(mat)
      }
      rownames(mat) <- var_names
      colnames(mat) <- obs_names
    }
    return(mat)
  }
}


# h5 to annotation
h5_to_df <- function(anno_gp_name){
  # recover the data frame
  df <- anno_gp_name[]
  rownames(df) <- df[['index']]
  df <- df[!names(df) == 'index']
  # recover the category or factor levels
  if(length(h5attr_names(anno_gp_name)>0)){
    for(k in h5attr_names(anno_gp_name)){
      cate_levels <- h5attr(anno_gp_name, k)
      k_cate <- cate_levels[df[[k]] + 1] # for begin 1
      df[[k]] <- factor(k_cate, levels = cate_levels)
    }
  }
  return(df)
}


# 
h5_to_adata <- function(h5 = NULL){
  if('annotation' %in% names(h5)){
    anno <- h5[['annotation']]
    if('observes' %in% names(anno)){
      obs_df <- h5_to_df(anno_gp_name = anno[['observes']])
      obsm <- rownames(obs_df)
      m1 <- 1
    }
    if('sdVariables' %in% names(anno)){
      sdvar_df <- h5_to_df(anno_gp_name = anno[['sdVariables']])
      sdvarm <- rownames(sdvar_df)
      m2 <- 2
    }
    if('ndVariables' %in% names(anno)){
      ndvar_df <- h5_to_df(anno_gp_name = anno[['ndVariables']])
      ndvarm <- rownames(ndvar_df)
      m3 <-3
    }
  }
  if(m1+m2+m3 == 6){
    if('normData' %in% names(h5) & !'rawData' %in% names(h5)){
      norm_data <- h5_to_matrix(gp_name = h5[['normData']], obs_names = obsm, var_names = ndvarm)
      seurat <- Seurat::CreateSeuratObject(counts = norm_data)
      # observes variable hvg
      seurat@meta.data <- obs_df
      seurat@assays$RNA@meta.features <- ndvar_df
      seurat@assays$RNA@var.features <- sdvarm
    } else if('rawData' %in% names(h5) & !'normData' %in% names(h5)){
      raw_data <- h5_to_matrix(gp_name = h5[['rawData']], obs_names = obsm, var_names = ndvarm)
      seurat <- Seurat::CreateSeuratObject(counts = raw_data)
      # observes variable hvg
      seurat@meta.data <- obs_df
      seurat@assays$RNA@meta.features <- ndvar_df
      seurat@assays$RNA@var.features <- sdvarm
    } else if('rawData' %in% names(h5) & 'normData' %in% names(h5)){
      norm_data <- h5_to_matrix(gp_name = h5[['normData']], obs_names = obsm, var_names = ndvarm)
      seurat <- Seurat::CreateSeuratObject(counts = norm_data)
      # observes variable hvg
      seurat@meta.data <- obs_df
      seurat@assays$RNA@meta.features <- ndvar_df
      seurat@assays$RNA@var.features <- sdvarm
      # raw data
      raw_data <- h5_to_matrix(gp_name = h5[['rawData']], obs_names = obsm, var_names = ndvarm)
      seurat@assays$RNA@counts <- raw_data
    }
    if('scaleData' %in% names(h5)){
      scale_data <- h5_to_matrix(gp_name = h5[['scaleData']], obs_names = obsm, var_names = sdvarm)
      if(class(scale_data) == 'matrix'){
        seurat@assays$RNA@scale.data <- scale_data
      }
    }
    if('dimReduction' %in% names(h5)){
      dimR <- h5[['dimReduction']]
      for(D in names(dimR)){
        d = tolower(D)
        dim_recover = t(dimR[[D]][,])
        dim_recover_ <- Seurat::CreateDimReducObject(embeddings = dim_recover, key = paste0(D, "_"), assay = "RNA")
        seurat@reductions[[d]] <- dim_recover_
      }
    }
    if('graphs' %in% names(h5)){
      graphs <- h5[["graphs"]]
      gra_list <- list(knn = "RNA_nn", snn = "RNA_snn")
      for(g in names(graphs)){
        graphs_neig_data = h5_to_matrix(gp_name = graphs[[g]], obs_names = obsm, var_names = obsm)
        seurat@graphs[[gra_list[[g]]]] <- graphs_neig_data
      }
    }
    if('metadata' %in% names(h5)){
      meta <- h5[['metadata/colors']]
      for( k in names(meta)){
        colo <- meta[[k]][]
        seurat@misc[[k]] <- colo
      }
    }
    return(seurat)
  }else{
    stop("Dataset lacks the necessary information")
  }
}

# read h5    
R_read_hdf5 <- function(file = NULL){
  if(is.null(file) | !file.exists(file)){
    stop('No such file or directory')
  } 
  h5 <- H5File$new(filename = file, mode = 'r')
  tryCatch({
    seurat <- h5_to_adata(h5 = h5 )
  }, error = function(err) {
    cat("error!", err, "\n")
  }, finally = {
    h5$close_all()
  }
  )
  return(seurat)
}