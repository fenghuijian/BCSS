

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

# write the hdf5
R_write_hdf5 <- function(data = NULL, 
                         file = NULL, 
                         unmodified.dataset = c("rawData", "normData", "scaleData", 
                                                "observes",  "variables", "PCA", "ICA", 
                                                "TSNE", "UMAP", "knn", "snn")
){
  if(is.null(data)){
    stop("object ", substitute(data), " not found")
  }
  else if(class(data) != "Seurat"){
    stop("oject ", substitute(data), " class is not Seurat object")
  }
  if(is.null(file)){
    stop("No such file or directory") 
  }
  h5 <- H5File$new(filename = file, mode = "a")
  # data is not null, data is Seurat object
  # file is not null, mode have to be unique
  if("dimNames" %in% names(h5)){
    if(!"dimNames" %in% unmodified.dataset){
      h5$link_delete("dimNames")
      dimNames <- h5$create_group("dimNames")
      dimNames[["obs_names"]] <- colnames(data@assays$RNA@data)
      dimNames[["var_names"]] <- rownames(data@assays$RNA@data)
      print("The old dimNames has been coveraged")
    }else{print("The dimNames has not been changed")}
  }
  else{
    dimNames <- h5$create_group("dimNames")
    dimNames[["obs_names"]] <- colnames(data@assays$RNA@data)
    dimNames[["var_names"]] <- rownames(data@assays$RNA@data)
    print("The new dimNames has been added")
  }
  # rewrite in rawData
  if("rawData" %in% names(h5)){
    if(!"rawData" %in% unmodified.dataset){
      h5$link_delete("rawData")
      rawData <- h5$create_group("rawData")
      rawData[["values"]] <- data@assays$RNA@counts@x
      rawData[["indices"]] <- data@assays$RNA@counts@i
      rawData[["indptr"]] <- data@assays$RNA@counts@p
      rawData[["dims"]] <- rev(data@assays$RNA@counts@Dim)
      print("The old rawData has been coveraged")
    }else{print("The rawData has not been changed")}
  } 
  else{
    # write in rawData
    rawData <- h5$create_group("rawData")
    rawData[["values"]] <- data@assays$RNA@counts@x
    rawData[["indices"]] <- data@assays$RNA@counts@i
    rawData[["indptr"]] <- data@assays$RNA@counts@p
    rawData[["dims"]] <- rev(data@assays$RNA@counts@Dim)
    print("The new rawData has been added")
  }
  # rewrite in normData
  if("normData" %in% names(h5)){
    if(!"normData" %in% unmodified.dataset){
      h5$link_delete("normData")
      normData <- h5$create_group("normData")
      normData[["values"]] <- data@assays$RNA@data@x
      normData[["indices"]] <- data@assays$RNA@data@i
      normData[["indptr"]] <- data@assays$RNA@data@p
      normData[["dims"]] <- rev(data@assays$RNA@data@Dim)
      print("The old normData has been coveraged")
    }else{print("The norm Data has not been changed")}
  }
  else{
    # write in normData
    normData <- h5$create_group("normData")
    normData[["values"]] <- data@assays$RNA@data@x
    normData[["indices"]] <- data@assays$RNA@data@i
    normData[["indptr"]] <- data@assays$RNA@data@p
    normData[["dims"]] <- rev(data@assays$RNA@data@Dim)
    print("The new normData has been added")
  }
  # rewrite in scaleData
  if("scaleData" %in% names(h5)){
    if(!"scaleData" %in% unmodified.dataset){
      h5$link_delete("scaleData")
      scaleData <- h5$create_group("scaleData")
      scaleData[["matrix"]] <- data@assays$RNA@scale.data
      scaleData[["dims"]] <- rev(dim(data@assays$RNA@scale.data))
      print("The old scaleData has been coveraged")
    }else{print("The scaleData has not been changed")}
  }
  else{
    # write in scaleData
    scaleData <- h5$create_group("scaleData")
    scaleData[["matrix"]] <- data@assays$RNA@scale.data
    scaleData[["dims"]] <- rev(dim(data@assays$RNA@scale.data))
    print("The new scaleData has been added")
  }
  # annotation
  if("annotation" %in% names(h5)){
    annotation <- h5[["annotation"]]
    if(!"observes" %in% unmodified.dataset){
      if(!"observes" %in% names(annotation) & !is.null(data@meta.data)){
        annotation[["observes"]] <- data@meta.data
        h5attr(annotation[["observes"]], "colnames") <- colnames(data@meta.data)
        print("The new annotation/observes has been added")
      }else{
        annotation$link_delete("observes")
        annotation[["observes"]] <- data@meta.data
        h5attr(annotation[["observes"]], "colnames") <- colnames(data@meta.data)
        print("The old annotation/observes has been coveraged")
      }
    }
    else{print("The annotation/observes has not been changed")}
    if(!"variables" %in% unmodified.dataset){
      if(!"variables" %in% names(annotation) & !is.null(data@assays$RNA@meta.features)){
        annotation[["variables"]] <- data@assays$RNA@meta.features
        h5attr(annotation[["variables"]], "colnames") <- colnames(data@assays$RNA@meta.features)
        print("The new annotation/variables has been added")
      }else{
        annotation$link_delete("variables")
        annotation[["variables"]] <- data@assays$RNA@meta.features
        h5attr(annotation[["variables"]], "colnames") <- colnames(data@assays$RNA@meta.features)
        print("The old annotation/variables has been coveraged")
      }
    }
    else{print("The annotation/variables has not been changed")}
  }else{
    annotation <- h5$create_group("annotation")
    # write in observes
    annotation[["observes"]] <- data@meta.data
    h5attr(annotation[["observes"]], "colnames") <- colnames(data@meta.data)
    print("The new annotation/observes has been added")
    # write in variables
    annotation[["variables"]] <- data@assays$RNA@meta.features
    h5attr(annotation[["variables"]], "colnames") <- colnames(data@assays$RNA@meta.features)
    print("The new annotation/variables has been added")
  }
  # dimReduction
  if("dimReduction" %in% names(h5)){
    dimReduction <- h5[["dimReduction"]]
    if(!"PCA" %in% unmodified.dataset){
      # rewrite in pca
      if(!"PCA" %in% names(dimReduction) & "pca" %in% names(data@reductions)){
        dimReduction[["PCA"]] <- data@reductions$pca@cell.embeddings
        print("The dimReduction/PCA has been added")
      }else if ("PCA" %in% names(dimReduction) & "pca" %in% names(data@reductions)){
        dimReduction$link_delete("PCA")
        dimReduction[["PCA"]] <- data@reductions$pca@cell.embeddings
        print("The old dimReduction/PCA has been coveraged")
      }
    }else{print("The dimReduction/PCA has not been changed")}
    if(!"ICA" %in% unmodified.dataset){
      # rewrite in ica
      if(!"ICA" %in% names(dimReduction) & "ica" %in% names(data@reductions)){
        dimReduction[["ICA"]] <- data@reductions$ica@cell.embeddings
        print("The dimReduction/ICA has been added")
      }else if("ICA" %in% names(dimReduction) & "ica" %in% names(data@reductions)){
        dimReduction$link_delete("ICA")
        dimReduction[["ICA"]] <- data@reductions$ica@cell.embeddings
        print("The old dimReduction/ICA has been coveraged")
      }
    }else{print("The dimReduction/ICA has not been changed")}
    if(!"TSNE" %in% unmodified.dataset){
      # rewrite in tsne
      if(!"TSNE" %in% names(dimReduction) & "tsne" %in% names(data@reductions)){
        dimReduction[["TSNE"]] <- data@reductions$tsne@cell.embeddings
        print("The dimReduction/TSNE has been added")
      }else if("TSNE" %in% names(dimReduction) & "tsne" %in% names(data@reductions)){
        dimReduction$link_delete("TSNE")
        dimReduction[["TSNE"]] <- data@reductions$tsne@cell.embeddings
        print("The old dimReduction/TSNE has been coveraged")
      }
    }else{print("The dimReduction/TSNE has not been changed")}
    if(!"UMAP" %in% unmodified.dataset){
      # rewrite in umap
      if(!"UMAP" %in% names(dimReduction) & "umap" %in% names(data@reductions)){
        dimReduction[["UMAP"]] <- data@reductions$umap@cell.embeddings
        print("The dimReduction/UMAP has been added")
      } else if("UMAP" %in% names(dimReduction) & "umap" %in% names(data@reductions)){
        dimReduction$link_delete("UMAP")
        dimReduction[["UMAP"]] <- data@reductions$umap@cell.embeddings
        print("The old dimReduction/UMAP has been coveraged")
      }
    }else{print("The dimReduction/UMAP has not been changed")}
  }else{
    dimReduction <- h5$create_group("dimReduction")
    if("pca" %in% names(data@reductions)){
      # write in pca
      dimReduction[["PCA"]] <- data@reductions$pca@cell.embeddings
      print("The new dimReduction/PCA has been added")
    }
    if("ica" %in% names(data@reductions)){
      # write in ica
      dimReduction[["ICA"]] <- data@reductions$ica@cell.embeddings
      print("The new dimReduction/ICA has been added")
    }
    if("tsne" %in% names(data@reductions)){
      # write in tsne
      dimReduction[["TSNE"]] <- data@reductions$tsne@cell.embeddings
      print("The new dimReduction/TSNE has been added")
    }
    if("umap" %in% names(data@reductions)){
      # write in umap
      dimReduction[["UMAP"]] <- data@reductions$umap@cell.embeddings
      print("The new dimReduction/UMAP has been added")
    }
  }
  if("graphs" %in% names(h5)){
    graphs <- h5[["graphs"]]
    if(!"knn" %in% unmodified.dataset){
      # rewrite in knn
      graphs$link_delete("knn")
      knn <- graphs$create_group("knn")
      knn[["values"]] <- data@graphs$RNA_nn@x
      knn[["indices"]] <- data@graphs$RNA_nn@i
      knn[["indptr"]] <- data@graphs$RNA_nn@p
      knn[["dims"]] <- data@graphs$RNA_nn@Dim
      print("The old graphs/knn has been coveraged")
    }else{print("The graphs/knn has not been changed")}
    if(!"snn" %in% unmodified.dataset){
      # rewrite in snn
      graphs$link_delete("snn")
      snn <- graphs$create_group("snn")
      snn[["values"]] <- data@graphs$RNA_snn@x
      snn[["indices"]] <- data@graphs$RNA_snn@i
      snn[["indptr"]] <- data@graphs$RNA_snn@p
      snn[["dims"]]  <- data@graphs$RNA_snn@Dim
      print("The old graphs/snn has been coveraged")
    }else{print("The graphs/snn has not been changed")}
  }else{
    graphs <- h5$create_group("graphs")
    if("RNA_nn" %in% names(data@graphs)){
      # write in knn
      knn <- graphs$create_group("knn")
      knn[["values"]] <- data@graphs$RNA_nn@x
      knn[["indices"]] <- data@graphs$RNA_nn@i
      knn[["indptr"]] <- data@graphs$RNA_nn@p
      knn[["dims"]] <- data@graphs$RNA_nn@Dim
      print("The new graphs/knn has been added")
    }
    if("RNA_snn" %in% names(data@graphs)){
      # rewrite in snn
      snn <- graphs$create_group("snn")
      snn[["values"]] <- data@graphs$RNA_snn@x
      snn[["indices"]] <- data@graphs$RNA_snn@i
      snn[["indptr"]] <- data@graphs$RNA_snn@p
      snn[["dims"]] <- data@graphs$RNA_snn@Dim
      print("The new graphs/snn has been added")
    }
  }
  # metaData
  h5$close_all()
}





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read the Hdf5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

R_read_hdf5 <- function(file = NULL){
  if(is.null(file) | !file.exists(file)){
    stop("No such file or directory")
  } 
  else if(file.exists(file)){
    h5 <- H5File$new(filename = file, mode = "r")
  }
  # obs_names and var names
  dimNames <- h5[["dimNames"]]
  obs_names <- dimNames[["obs_names"]][]
  var_names <- dimNames[["var_names"]][]
  # raw Data recover
  rawData <- h5[["rawData"]]
  raw_data <- sparseMatrix(i = rawData[["indices"]][], 
                           p = rawData[["indptr"]][], 
                           x = rawData[["values"]][],
                           dims = rev(rawData[["dims"]][]),
                           dimnames = list(var_names, obs_names), 
                           index1 = FALSE)
  seurat <- CreateSeuratObject(counts = raw_data)
  # norm Data recover
  normData <- h5[["normData"]]
  norm_data <- sparseMatrix(i = normData[["indices"]][], 
                            p = normData[["indptr"]][], 
                            x = normData[["values"]][],
                            dims = rev(normData[["dims"]][]), 
                            dimnames = list(var_names, obs_names), 
                            index1 = FALSE)
  seurat@assays$RNA@data <- norm_data
  # scale Data recover
  scaleData <- h5[["scaleData"]]
  scale_data <- scaleData[["matrix"]][,]
  rownames(scale_data) <- var_names
  colnames(scale_data) <- obs_names
  seurat@assays$RNA@scale.data <- scale_data
  # meta data recover 
  annotation <- h5[["annotation"]]
  # 
  metadata <- annotation[["observes"]][]
  rownames(metadata) <- obs_name
  seurat@meta.data <- metadata
  features <- annotation[["variables"]][]
  seurat@assays$RNA@meta.features <- features
  # dimReduction
  #pca
  dimReduction <- h5[["dimReduction"]]
  if("PCA" %in% names(dimReduction)){
    pca <- dimReduction[["PCA"]][,]
    rownames(pca) <- obs_names
    pca_ <- CreateDimReducObject(embeddings = pca, key = "PCA_")
    seurat@reductions$pca <- pca_
  }
  if("ICA" %in% names(dimReduction)){
    ica <- dimReduction[["ICA"]][,]
    rownames(ica) <- obs_names
    ica_ <- CreateDimReducObject(embeddings = ica, key = "ICA_")
    seurat@reductions$ica <- ica_
  }
  if("TSNE" %in% names(dimReduction)){
    tsne <- dimReduction[["TSNE"]][,]
    rownames(tsne) <- obs_name
    tsne_ <- CreateDimReducObject(embeddings = tsne, key = "TSNE_")
    seurat@reductions$tsne <- tsne_
  }
  if("UMAP" %in% names(dimReduction)){
    umap <- dimReduction[["UMAP"]][,]
    rownames(umap) <- obs_names
    umap_ <- CreateDimReducObject(embeddings = umap, key = "UMAP_")
    seurat@reductions$umap <- umap_
  }
  # graph
  graphs <- h5[["graphs"]]
  if("knn" %in% names(graphs)){
    knn <- sparseMatrix(i = graphs[["knn/indices"]][], 
                        p = graphs[["knn/indptr"]][], 
                        x = graphs[["knn/values"]][], 
                        dims = graphs[["knn/dims"]][], 
                        dimnames = list(obs_names, obs_names), 
                        index1 = FALSE)
    seurat@graphs$RNA_nn <- knn
  }
  if("snn" %in% names(graphs)){
    snn <- sparseMatrix(i = graphs[["snn/indices"]][], 
                        p = graphs[["snn/indptr"]][], 
                        x = graphs[["snn/values"]][], 
                        dims = graphs[["snn/dims"]][], 
                        dimnames = list(obs_names, obs_names),
                        index1 = FALSE)
    seurat@graphs$RNA_snn <- snn
  }
  # 
  return(seurat)
}
























