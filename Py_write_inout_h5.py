# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 23:23:03 2019

@author: fenghuijian
"""

import h5py
import numpy as np
import scipy
from scipy import sparse
import scanpy as sc
import pandas as pd
import anndata
from pandas.api.types import is_string_dtype, is_categorical



# dataframe to the array for the compound datatypes hdf5
def df_to_arrays(df):
    uns={}
    names = ["index"]
    if is_string_dtype(df.index):
        index = df.index.values.astype(h5py.special_dtype(vlen=str))
    else:
        index = df.index.values
    arrays = [index]
    for k in df.keys():
        names.append(k)
        if is_string_dtype(df[k]) and not is_categorical(df[k]):
            arrays.append(df[k].values.astype(h5py.special_dtype(vlen=str)))
        elif is_categorical(df[k]):
            arrays.append(df[k].cat.codes)
            uns[k + "_category"] = df[k].cat.categories
        else:
            arrays.append(df[k].values)
    dt = [d.dtype for d in arrays]
    return np.rec.fromarrays(arrays, dtype = {"names" : names, "formats" : dt}), uns




# write hdf5
def Py_write_hdf5(data=None, rawdata=None, save_basic = False, file = None):
    def namestr(obj, namespace):
        return [name for name in namespace if namespace[name] is obj]
    if file is None:
        raise OSError("No such file or directory")
    # w Create file, truncate if exists
    h5 = h5py.File(name=file, mode="w")   
    if not isinstance(data, anndata.AnnData):
        h5.close()
        raise TypeError("object '%s' class is not anndata.AnnData object" %namestr(data, globals())[0])
    else:
        # normData must be corresponded data.raw.X
        normData=h5.create_group("normData")
        ndata_i=normData.create_dataset("indices", data=data.raw.X.indices, dtype="i4")
        ndata_p=normData.create_dataset("indptr", data=data.raw.X.indptr, dtype="i4")
        ndata_x=normData.create_dataset("values", data=data.raw.X.data, dtype="f4")
        ndata_dims=normData.create_dataset("dims", data=data.raw.X.shape)
        # scaleData must be corrsesponded data.X
        scaleData=h5.create_group('scaleData')
        sdata_matrix=scaleData.create_dataset("matrix", data=data.X, dtype="f4")
        sdata_dims=scaleData.create_dataset("dims", data=data.X.shape) 
        # annotation
        annotation=h5.create_group('annotation')
        # -- observes
        obs,obs_uns = df_to_arrays(df = data.obs)
        observes = annotation.create_dataset("observes", data = obs)
        if not obs_uns is None:
            for o in obs_uns:
                observes.attrs[o]=obs_uns[o].values.astype(h5py.special_dtype(vlen=str))
        # -- variables
        var,var_uns = df_to_arrays(df = data.var)
        variables = annotation.create_dataset("variables", data = var)
        if not var_uns is None:
            for v in var_uns:
                variables.attrs[v]=var_uns[v].values.astype(h5py.special_dtype(vlen=str))
        # data.var.to_hdf(file, key = "annotation/variables", mode = "w")
        if save_basic is True:
            if isinstance(rawdata, anndata.AnnData):
                # rawData must be corresponded rawdata.X
                rawData=h5.create_group("rawData")
                rdata_i=rawData.create_dataset("indices", data=rawdata.X.indices, dtype="i4")
                rdata_p=rawData.create_dataset("indptr", data=rawdata.X.indptr, dtype="i4")
                rdata_x=rawData.create_dataset("values", data=rawdata.X.data, dtype="f4")
                rdata_dims=rawData.create_dataset("dims", data=rawdata.X.shape)
            else:
                h5.close()
                raise TypeError("object '%s' class is not anndata.AnnData object" %namestr(rawdata, globals())[0])
        else:
            # {"X_pca":"PCA", "X_ica":"ICA", "X_nmf":"NMF", "X_tsne":"TSNE", "X_umap":"UMAP", "X_da":"DC"}
            dimReduction=h5.create_group("dimReduction")
            dr={"X_pca":"PCA", "X_ica":"ICA", "X_nmf":"NMF", "X_tsne":"TSNE", "X_umap":"UMAP", "X_dc":"DC"}
            for k in [k for k in data.obsm.keys()]:
                    dimReduction.create_dataset(dr[k], data=data.obsm[k], dtype="f4")          
            # graphs
            graphs=h5.create_group("graphs")
            # -- knn like
            knn=graphs.create_group("knn")
            knn_i=knn.create_dataset("indices", data = data.uns['neighbors']['distances'].indices, dtype="i4")
            knn_p=knn.create_dataset("indptr", data = data.uns['neighbors']['distances'].indptr, dtype="i4")
            knn_x=knn.create_dataset("values", data = data.uns['neighbors']['distances'].data, dtype="f4")
            knn_i=knn.create_dataset("dims", data = data.uns['neighbors']['distances'].shape)
            # -- snn like
            snn=graphs.create_group("snn")
            snn_i=snn.create_dataset("indices", data = data.uns['neighbors']['connectivities'].indices, dtype="i4")
            snn_p=snn.create_dataset("indptr", data = data.uns['neighbors']['connectivities'].indptr, dtype="i4")
            snn_x=snn.create_dataset("values", data = data.uns['neighbors']['connectivities'].data, dtype="f4")
            snn_i=snn.create_dataset("dims", data = data.uns['neighbors']['connectivities'].shape)
            
            # close 
            h5.close()
    return "The data was written in successfully"
    


file = "py_test_2.hdf5"
file = "py_test.hdf5"
Py_write_hdf5(data = adata, file = file)
h5.close()

add=Py_read_hdf5(file)
############ read the hdf5 file
def Py_read_hdf5(file = None, use_raw = False):
    def h5_to_csr(gp):
        x = gp["values"][()]
        indices = gp["indices"][()]
        indptr = gp["indptr"][()]
        shape = gp["dims"][()]
        sp_csr = sparse.csr_matrix((x,indices,indptr),shape=shape, dtype = "f4")
        return sp_csr
    def h5_to_df(anno_gp):
        dict_ov={"observes":"obs_names", "variables":"var_names"}
        list=[]
        for k in dict_ov:
            dict_ce={}
            data = anno_gp[k][()]
            for j in data.dtype.names:
                if j=="index":
                    dict_ce[dict_ov[k]] = data[j]
                else:
                    dict_ce[j] = data[j]
            list.append(dict_ce)
        return list
    if file is None:
        raise OSError("No such file or directory")
    h5 = h5py.File(name=file, mode="r")
    if use_raw:
        rdata=h5["rawData"]
        annotation=h5["annotation"]
        sp_csr=h5_to_csr(rdata)
        anno=h5_to_df(annotation)
        adata = anndata.AnnData(sp_csr, anno[0], anno[1], dtype=sp_csr.dtype.name)
    else:
        ndata=h5["normData"]
        sdata=h5["scaleData/matrix"][()]
        dimR=h5["dimReduction"]
        graphs=h5["graphs"]
        annotation=h5["annotation"]
        sp_csr=h5_to_csr(ndata)
        anno=h5_to_df(annotation)
        adata = anndata.AnnData(sp_csr, anno[0], anno[1], dtype=sp_csr.dtype.name)
        adata.raw = adata
        adata.X = sdata
        dr={"PCA":"X_pca", "ICA":"X_ica", "NMF":"X_nmf", "TSNE": "X_tsne", "UMAP": "X_umap", "DC": "X_dc"}
        for k in dimR.keys():
            adata.obsm[dr[k]] = dimR[k][()]
        # graphs :
        neig = {"knn":"distances", "snn":"connectivities"}
        neig_dict ={}
        for g in graphs.keys():
            neig_dict[neig[g]]= h5_to_csr(graphs[g])
        adata.uns["neighbors"] = neig_dict
    return adata









