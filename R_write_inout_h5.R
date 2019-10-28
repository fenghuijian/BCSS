# library the prerequisite packages
library(Seurat)
library(hdf5r)
library(Matrix)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Write the Hdf5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# matrix to h5
matrix_to_h5 <- function(mat, h5_gp, gp_name = NULL){
  if(!gp_name %in% names(h5_gp)){
    h5mat = h5_gp$create_group(gp_name)
  }
  else{
    h5mat = h5_gp[[gp_name]]
  }
  if('dgCMatrix' %in% class(mat)){
    h5mat[['values']] <- slot(object = mat, name = 'x')
    h5mat[['indices']] <- slot(object = mat, name = 'i')
    h5mat[['indptr']] <- slot(object = mat, name = 'p')
    h5mat[['dims']] <- rev(slot(object = mat, name = 'Dim'))
    h5attr(h5mat, 'datatype') <- 'SparseMatrix'
  }
  else if('matrix' %in% class(mat)){
    h5mat[['matrix']] <- mat
    h5mat[['dims']] <- rev(dim(mat))
    h5attr(h5mat, 'datatype') <- 'Array'
  }
  else if('Graph' %in% class(mat)){
    h5mat[['values']] <- slot(object = mat, name = 'x')
    h5mat[['indices']] <- slot(object = mat, name = 'i')
    h5mat[['indptr']] <- slot(object = mat, name = 'p')
    h5mat[['dims']] <- rev(slot(object = mat, name = 'Dim'))
    h5attr(h5mat, 'datatype') <- 'SparseMatrix'
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

# seurat be converted to h5
adata_to_h5 <- function(adata=NULL, h5=NULL, assay.name = NULL){
  if(assay.name %in% names(slot(object = adata, name = 'assays'))){
    slot_assay <- slot(object = adata, name = 'assays')[[assay.name]]
    if(all(slot_assay@counts@x == slot_assay@data@x)){
      # only have rdata (rdata == ndata)
      rdata <- slot(object = slot_assay, name = 'counts')
      matrix_to_h5(mat = rdata, h5_gp = h5, gp_name = 'rawData')
      
    } else{
      # have rdata and ndata (rdata != ndata)
      rdata <- slot(object = slot_assay, name = 'counts')
      matrix_to_h5(mat = rdata, h5_gp = h5, gp_name = 'rawData')
      ndata <- slot(object = slot_assay, name = 'data')
      matrix_to_h5(mat = ndata, h5_gp = h5, gp_name = 'normData')
    }
    if(length(slot(slot_assay, 'scale.data'))>0){
      sdata <- slot(object = slot_assay, name = 'scale.data')
      matrix_to_h5(mat = sdata, h5_gp = h5, gp_name = 'scaleData')
    }else{
      sdata <- NULL
    }
    # anno
    annotation <- h5$create_group('annotation')
    # -- observes
    obs = slot(object = adata, name = 'meta.data')
    obs[['index']] <- rownames(obs)
    df_to_h5(df = obs, h5_anno = annotation, anno_dataset = 'observes')
    
    # -- variables
    ndvar <- slot(object = slot_assay, name = 'meta.features')
    hvg <- slot(object = slot_assay, name = 'var.features')
    ndvar[['index']] <- rownames(ndvar)
    ndvar[['highly_variable']] <- ndvar[['index']] %in% hvg
    df_to_h5(df = ndvar, h5_anno = annotation, anno_dataset = 'ndVariables')
    if(!is.null(sdata)){
      if(all(nrow(sdata) == nrow(ndata))){
        df_to_h5(df = ndvar, h5_anno = annotation, anno_dataset = 'sdVariables')
      } else{
        sdvarm <- rownames(sdata)
        sdvar <- ndvar[sdvarm, ]
        df_to_h5(df = sdvar, h5_anno = annotation, anno_dataset = 'sdVariables')
      }
    }
    if(length(adata@reductions)>0){
      dimReduction <- h5$create_group('dimReduction')
      dimr <- slot(object = adata, name = 'reductions')
      for(d in names(dimr)){
        D = toupper(d)
        dimReduction[[D]] <- t(slot(dimr[[d]],'cell.embeddings'))
      }
    }
    if(length(adata@graphs)>0){
      graph_df <- slot(object = adata, 'graphs')
      if(length(grep(assay.name, names(graph_df))) == 2){
        graphs <- h5$create_group('graphs')
        gra_list <- list()
        gra_list[[paste0(assay.name, '_nn')]] <- 'knn'
        gra_list[[paste0(assay.name, '_snn')]] <- 'snn'
        for(g in names(graph_df)){
          matrix_to_h5(mat=graph_df[[g]], h5_gp = graphs, gp_name = gra_list[[g]])
        }
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
  }else{
    stop('Please enter the correct assay name')
  }
  
}

# R write hdf5
R_write_hdf5 <- function(adata = NULL, file = NULL, assay.name = 'RNA'){
  if(is.null(file)){
    stop('No such file or directory') 
  }
  if(class(adata) != 'Seurat'){
    stop('oject ', substitute(adata), ' class is not Seurat object')
  }
  h5 <- H5File$new(filename = file, mode = 'w')
  tryCatch({
    adata_to_h5(adata = adata, h5 = h5, assay.name = assay.name)
    h5attr(h5, 'assay_name') <- assay.name
  }, error = function(e) {
    print(e)
  }, finally = {
    h5$close_all()
  }
  )
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read the Hdf5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# h5 to matrix
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





# h5 to the seurat
h5_to_adata <- function(h5 = NULL, assay.name = NULL){
  if(h5attr(h5, 'assay_name') == assay.name){
    if('annotation' %in% names(h5)){
      anno <- h5[['annotation']]
      if('observes' %in% names(anno)){
        obs_df <- h5_to_df(anno_gp_name = anno[['observes']])
        obsm <- rownames(obs_df)
      }else{
        obs_df <- NULL
      }
      if('sdVariables' %in% names(anno)){
        sdvar_df <- h5_to_df(anno_gp_name = anno[['sdVariables']])
        sdvarm <- rownames(sdvar_df)
      }else{
        sdvar_df <- NULL
      }
      if('ndVariables' %in% names(anno)){
        ndvar_df <- h5_to_df(anno_gp_name = anno[['ndVariables']])
        ndvarm <- rownames(ndvar_df)
      }else{
        ndvar_df <- NULL
      }
    }
    if('normData' %in% names(h5) & !'rawData' %in% names(h5)){
      # only have norm data
      if(!is.null(obs_df) & !is.null(ndvar_df)){
        norm_data <- h5_to_matrix(gp_name = h5[['normData']], obs_names = obsm, var_names = ndvarm)
        seurat <- Seurat::CreateSeuratObject(counts = norm_data, assay = assay.name)
        slot(object = seurat, name = 'meta.data') <- obs_df
        seurat@assays[[assay.name]]@meta.features <- ndvar_df
      } else{stop('Data structures are lack of cellular or genetic annotation')}
    } else if('rawData' %in% names(h5) & !'normData' %in% names(h5)){
      # only have raw counts
      if(!is.null(obs_df) & !is.null(ndvar_df)){
        raw_data <- h5_to_matrix(gp_name = h5[['rawData']], obs_names = obsm, var_names = ndvarm)
        seurat <- Seurat::CreateSeuratObject(counts = raw_data, assay = assay.name)
        slot(object = seurat, name = 'meta.data') <- obs_df
        seurat@assays[[assay.name]]@meta.features <- ndvar_df
      } else{stop('Data structures are lack of cellular or genetic annotation')}
    } else if('rawData' %in% names(h5) & 'normData' %in% names(h5)){
      if(!is.null(obs_df) & !is.null(ndvar_df)){
        norm_data <- h5_to_matrix(gp_name = h5[['normData']], obs_names = obsm, var_names = ndvarm)
        # raw data
        raw_data <- h5_to_matrix(gp_name = h5[['rawData']], obs_names = obsm, var_names = ndvarm)
        seurat <- Seurat::CreateSeuratObject(counts = raw_data, assay = assay.name)
        slot(object = seurat, name = 'meta.data') <- obs_df
        seurat@assays[[assay.name]]@data <- norm_data
        seurat@assays[[assay.name]]@meta.features <- ndvar_df
      } else{stop('Data structures are lack of cellular or genetic annotation')}
    } else{stop('Problem with data structure')}
    if('scaleData' %in% names(h5)){
      if(!is.null(sdvar_df)){
        scale_data <- h5_to_matrix(gp_name = h5[['scaleData']], obs_names = obsm, var_names = sdvarm)
        seurat@assays[[assay.name]]@scale.data <- scale_data
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
  }else{stop('The assay information of h5 is inconsistent')}
  
}

# Read h5   
R_read_hdf5 <- function(file = NULL, assay.name = 'RNA'){
  if(is.null(file) | !file.exists(file)){
    stop('No such file or directory')
  } 
  h5 <- H5File$new(filename = file, mode = 'r')
  tryCatch({
    seurat <- h5_to_adata(h5 = h5, assay.name = assay.name)
  }, error = function(e) {
    print(e)
  }, finally = {
    h5$close_all()
  }
  )
  return(seurat)
}








