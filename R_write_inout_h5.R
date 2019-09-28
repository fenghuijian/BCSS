

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
h5 <- H5File$new(file, "w")

file ="pd.h5"

setwd("B:/Provincial_laboratory")
h5$close_all()
# write the hdf5
R_write_hdf5 <- function(data = NULL, 
                         file = NULL
){
  if(is.null(file)){
    stop("No such file or directory") 
  }
  h5 <- H5File$new(filename = file, mode = "w")
  if(class(data) != "Seurat"){
    h5$close_all()
    stop("oject ", substitute(data), " class is not Seurat object")
  }
  # data is not null, data is Seurat object
  # file is not null, mode have to be unique
  # write in rawData
  rawData <- h5$create_group("rawData")
  rdata <- GetAssayData(object = data, slot = "data", assay = "RNA")
  rawData[["values"]] <- rdata@x
  rawData[["indices"]] <- rdata@i
  rawData[["indptr"]] <- rdata@p
  rawData[["dims"]] <- rev(rdata@Dim)
  # write in normData
  normData <- h5$create_group("normData")
  ndata <- GetAssayData(object = data, slot = "data", assay = "RNA")
  normData[["values"]] <- ndata@x
  normData[["indices"]] <- ndata@i
  normData[["indptr"]] <- ndata@p
  normData[["dims"]] <- rev(ndata@Dim)
  # write in scaleData
  if(!is.null(data@assays$RNA@scale.data)){
    scaleData <- h5$create_group("scaleData")
    sdata <- GetAssayData(object = data, slot = "scale.data", assay = "RNA")
    scaleData[["matrix"]] <- sdata
    scaleData[["dims"]] <- rev(dim(sdata))
  }
  # annotation
  annotation <- h5$create_group("annotation")
  # write in observes
  anno_obs <- data@meta.data
  anno_obs$index <- rownames(anno_obs)
  annotation[["observes"]] <- anno_obs
  anno_list <- list(louvain_category = "seurat_clusters")
  for(o in names(anno_list)){
    h5attr(annotation[["observes"]], o) <- levels(data@meta.data[, anno_list[[o]]])
  }
  # write in variables
  anno_var <- data@assays$RNA@meta.features
  anno_var$gene_ids <- rownames(anno_var)
  annotation[["variables"]] <- anno_var
  # dimReduction
  dimReduction <- h5$create_group("dimReduction")
  dim_list <- list(pca = "PCA", ica = "ICA", nmf = "NMF", tsne = "TSNE", 
                   umap = "UMAP", dc = "DC")
  for(d in names(data@reductions)){
    dimReduction[[dim_list[[d]]]] <- t(data@reductions[[d]]@cell.embeddings)
  }
  # graphs
  graphs <- h5$create_group("graphs")
  gra_list <- list(RNA_nn = "knn", RNA_snn = "snn")
  for(g in names(data@graphs)){
    ga <- graphs$create_group(gra_list[[g]])
    ga[["values"]] <- data@graphs[[g]]@x
    ga[["indices"]] <- data@graphs[[g]]@i
    ga[["indptr"]] <- data@graphs[[g]]@p
    ga[["dims"]] <- rev(data@graphs[[g]]@Dim)
  }
  # metaData
  h5$close_all()
  return("The data was written in successfully")
}

R_write_hdf5(data = pbmc, file = file)



h5 <- H5File$new(file, "r")

head(h5[["annotation/observes"]][])


h5$close_all()

a = pbmc@assays$RNA@scale.data


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

h5$ls()
h7 <- H5File$new("C:/Users/fenghuijian/py_test_2.hdf5", "r")
h6 <- H5File$new("B:/Provincial_laboratory/pd.h5", "r")
h6$close_all()


annotation$create_dataset(n)

h6[["annotation/observes"]]
h7[["annotation/observes"]]





h5 <- H5File$new("pdd.h5", "w")

annotation=h5$create_group("annotation")
logical_example <- H5T_LOGICAL$new(include_NA = TRUE)
## we could also use h5types$H5T_LOGICAL or h5types$H5T_LOGICAL_NA
logical_example$get_labels()
logical_example$get_values()
cpd_example <- H5T_COMPOUND$new(c("Double_col", "Int_col", "Logical_col"), 
                                dtypes = list(h5types$H5T_NATIVE_DOUBLE, 
                                              h5types$H5T_NATIVE_INT, logical_example))
cpd_example

test <- annotation$create_dataset(name = "test", dtype = cpd_example, space = )

uint2_dt <- h5types$H5T_NATIVE_UINT32$set_size(1)$set_precision(2)$set_sign(h5const$H5T_SGN_NONE)
ds_create_pl_nbit <- H5P_DATASET_CREATE$new()
ds_create_pl_nbit$set_chunk(c(10, 10))$set_fill_value(uint2_dt, 1)$set_nbit()

c("da", 1)

