

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
                                                "TSNE", "UMAP", "knn", "snn"),
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
  
  # rewrite in rawData
  if("rawData" %in% names(h5)){
    if(!"rawData" %in% unmodified.dataset){
      rawData <- h5[["rawData"]]
      rawData[["values"]][] <- data@assays$RNA@counts@x
      rawData[["indices"]][] <- data@assays$RNA@counts@i
      rawData[["indptr"]][] <- data@assays$RNA@counts@p
      rawData[["dims"]][] <- rev(data@assays$RNA@counts@Dim)
      print("The rawData has been changed")
    }
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
      normData <- h5[["normData"]]
      normData[["values"]][] <- data@assays$RNA@data@x
      normData[["indices"]][] <- data@assays$RNA@data@i
      normData[["indptr"]][] <- data@assays$RNA@data@p
      normData[["dims"]][] <- rev(data@assays$RNA@data@Dim)
      print("The normData has been changed")
    }
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
      scaleData <- h5[["scaleData"]]
      scaleData[["matrix"]][,] <- data@assays$RNA@scale.data
      scaleData[["dims"]][] <- rev(dim(data@assays$RNA@scale.data))
      print("The scaleData has been changed")
    }
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
      # rewrite in observes
      annotation[["observes"]][] <- data@meta.data
      h5attr(annotation[["observes"]], "rownames") <- rownames(data@meta.data)
      h5attr(annotation[["observes"]], "colnames") <- colnames(data@meta.data)
      print("The annotation/observes has been changed")
    }
    if(!"variables" %in% unmodified.dataset){
      # rewrite in observes 
      annotation[["variables"]][] <- data.frame(gene_id = colnames(data@assays$RNA@data))
      h5attr(annotation[["variables"]], "rownames") <- colnames(data@assays$RNA@data)
      print("The annotation/variables has been changed")
    }
  }else{
    annotation <- h5$create_group("annotation")
      # write in observes
    annotation[["observes"]] <- data@meta.data
    h5attr(annotation[["observes"]], "rownames") <- rownames(data@meta.data)
    h5attr(annotation[["observes"]], "colnames") <- colnames(data@meta.data)
    print("The new annotation/observes has been added")
      # write in variables
    annotation[["variables"]] <- data.frame(gene_id = colnames(data@assays$RNA@data))
    h5attr(annotation[["variables"]], "rownames") <- colnames(data@assays$RNA@data)
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
        dimReduction[["PCA"]][,] <- data@reductions$pca@cell.embeddings
        print("The dimReduction/PCA has been changed")
      }
    }
    if(!"ICA" %in% unmodified.dataset){
      # rewrite in ica
      if(!"ICA" %in% names(dimReduction) & "ica" %in% names(data@reductions)){
        dimReduction[["ICA"]] <- data@reductions$ica@cell.embeddings
        print("The dimReduction/ICA has been added")
      }else if("ICA" %in% names(dimReduction) & "ica" %in% names(data@reductions)){
        dimReduction[["ICA"]][,] <- data@reductions$ica@cell.embeddings
        print("The dimReduction/ICA has been changed")
      }
    }
    if(!"TSNE" %in% unmodified.dataset){
      # rewrite in tsne
      if(!"TSNE" %in% names(dimReduction) & "tsne" %in% names(data@reductions)){
        dimReduction[["TSNE"]] <- data@reductions$tsne@cell.embeddings
        print("The dimReduction/TSNE has been added")
      }else if("TSNE" %in% names(dimReduction) & "tsne" %in% names(data@reductions)){
        dimReduction[["TSNE"]][,] <- data@reductions$tsne@cell.embeddings
        print("The dimReduction/TSNE has been changes")
      }
    }
    if(!"UMAP" %in% unmodified.dataset){
      # rewrite in umap
      if(!"UMAP" %in% names(dimReduction) & "umap" %in% names(data@reductions)){
        dimReduction[["UMAP"]] <- data@reductions$umap@cell.embeddings
        print("The dimReduction/UMAP has been added")
      } else if("UMAP" %in% names(dimReduction) & "tsne" %in% names(data@reductions)){
        dimReduction[["UMAP"]][,] <- data@reductions$umap@cell.embeddings
        print("The dimReduction/UMAP has been changed")
      }
    }
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
      knn <- graphs[["knn"]]
      knn[["values"]][] <- data@graphs$RNA_nn@x
      knn[["indices"]][] <- data@graphs$RNA_nn@i
      knn[["indptr"]][] <- data@graphs$RNA_nn@p
      knn[["dims"]][] <- data@graphs$RNA_nn@Dim
      print("The graphs/knn has been changed")
    }
    if(!"snn" %in% unmodified.dataset){
      # rewrite in snn
      snn <- graphs[["snn"]]
      snn[["values"]][] <- data@graphs$RNA_snn@x
      snn[["indices"]][] <- data@graphs$RNA_snn@i
      snn[["indptr"]][] <- data@graphs$RNA_snn@p
      snn[["dims"]][] <- data@graphs$RNA_snn@Dim
      print("The graphs/snn has been changed")
    }
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