# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 23:23:03 2019

@author: fenghuijian
"""
# import the ne
import h5py
import numpy as np
import scipy
from scipy import sparse
import scanpy as sc
import pandas as pd
import anndata
from pandas.api.types import is_string_dtype, is_categorical
#import re
#import os


# dataframe to the array for the compound datatypes hdf5
def df_to_h5(df, h5_anno, anno_dataset=None, anno_gp_name=None, anno_gp_dataset=None):
    #to array
    cate={}
    names = ['index']
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
            cate[k] = df[k].cat.categories
        else:
            arrays.append(df[k].values)
    dt = [d.dtype for d in arrays]
    h5_df = np.rec.fromarrays(arrays, dtype = {'names' : names, 'formats' : dt})
    # to h5
    if not anno_gp_name:
        h5_anno_ds = h5_anno.create_dataset(anno_dataset, data=h5_df)
        for o in cate:
            h5_anno_ds.attrs[o] = cate[o].values.astype(h5py.special_dtype(vlen=str))
    else:
        if anno_gp_name not in h5_anno.keys():
            h5_anno_gp = h5_anno.create_group(anno_gp_name)
        else: 
            h5_anno_gp = h5_anno[anno_gp_name]
        h5_anno_gp_ds = h5_anno_gp.create_dataset(anno_gp_dataset, data=h5_df)
        for p in cate:
            h5_anno_gp_ds.attrs[p] = cate[p].values.astype(h5py.special_dtype(vlen=str))
    return 


# glabol function
def namestr(obj, namespace):
    return [name for name in namespace if namespace[name] is obj]


# numpy array or scipy sparse matrix save to the h5
def matrix_to_h5(mat, h5_gp, gp_name = None):
    if gp_name not in h5_gp.keys():
        h5mat = h5_gp.create_group(gp_name)
    else:
        h5mat = h5_gp[gp_name]
    if isinstance(mat, scipy.sparse.csr.csr_matrix):
        h5mat_i = h5mat.create_dataset("indices", data=mat.indices)
        h5mat_p = h5mat.create_dataset("indptr", data=mat.indptr)
        h5mat_x = h5mat.create_dataset("values", data=mat.data, dtype=np.float32)
        h5mat_dims = h5mat.create_dataset("dims", data=mat.shape)
        h5mat.attrs["datatype"] = "SparseMatrix"
    elif isinstance(mat, np.ndarray):
        h5mat_mat = h5mat.create_dataset("matrix", data=mat, dtype=np.float32)
        h5mat_dims = h5mat.create_dataset("dims", data=mat.shape)
        h5mat.attrs['datatype'] = 'Array'
    return 




# write hdf5
def Py_write_hdf5(adata=None, file = None):
    if file is None:
        raise OSError("No such file or directory")
    if not isinstance(adata, anndata.AnnData):
        raise TypeError("object '%s' class is not anndata.AnnData object" %namestr(adata, globals())[0])
    # w Create file, truncate if exists
    h5 = h5py.File(name=file, mode="w")   
    try:
        # raw counts data is located the adata.layers['counts']
        if "counts" in adata.layers.keys():
            rdata = adata.layers["counts"]
            # write in
            matrix_to_h5(mat=rdata, h5_gp=h5, gp_name="rawData")
        # normalized data is located the adata.raw.X
        ndata = adata.raw.X
        matrix_to_h5(mat=ndata, h5_gp=h5, gp_name="normData")
        # scaledata is located the data.X
        sdata = adata.X
        matrix_to_h5(mat=sdata, h5_gp=h5, gp_name="scaleData")
        # annotation
        annotation = h5.create_group('annotation')
        df_to_h5(df=adata.obs, h5_anno=annotation, anno_dataset='observes')
        # -- variables # ---- scaledata - variables # ---- normdata - variables
        df_to_h5(df=adata.var, h5_anno=annotation, anno_gp_name='variables', anno_gp_dataset='sdVariables')
        df_to_h5(df=adata.raw.var, h5_anno=annotation, anno_gp_name='variables', anno_gp_dataset='ndVariables')
        # dimension reduction
        dimReduction=h5.create_group("dimReduction")
        for k in [k for k in adata.obsm.keys()]:
            K = re.sub("^.*_", "", k).upper()
            dimReduction.create_dataset(K, data=adata.obsm[k], dtype=np.float32)          
        # graphs
        gra_dict={"distances":"knn", "connectivities":"snn"}
        graph_data = adata.uns['neighbors']
        # -- save the neighbor graphs
        graphs = h5.create_group('graphs')
        for g in gra_dict.keys():
            matrix_to_h5(mat=graph_data[g], h5_gp=graphs, gp_name=gra_dict[g])
    except Exception as e:
        print('Error:',e)
    finally:
        h5.close()
    return 





############ read the hdf5 file

# h5 to the scipy sparse csr format and ndarray
def h5_to_matrix(gp_name):
    if gp_name.attrs['datatype'] == 'SparseMatrix':
        x = gp_name["values"][()].astype(np.float32)
        indices = gp_name["indices"][()]
        indptr = gp_name["indptr"][()]
        shapes = gp_name["dims"][()]
        mat = sparse.csr_matrix((x,indices,indptr),shape=shapes, dtype = np.float32)
    elif gp_name.attrs['datatype'] == 'Array':
        mat = gp_name['matrix'][()].astype(np.float32)
    return mat





# h5 to the pandas dataframe
def h5_to_df(anno_gp_name):
    # h5 to pandas dataframe
    dict_ce={}
    data = anno_gp_name[()]
    for j in data.dtype.names:
        if data[j].dtype == np.object:
            if j == "index":
                dict_ce["index"] = np.array(data[j].astype("str").tolist(), dtype =np.object0)
            else:
                dict_ce[j] = np.array(data[j].astype("str").tolist(), dtype =np.object0)
        else:
            dict_ce[j] = data[j]
    df = pd.DataFrame(dict_ce)
    df = df.set_index('index')
    # recover pandas the category
    if len(anno_gp_name.attrs.keys()) > 0:
        for k in anno_gp_name.attrs.keys():
            cate_levels = np.array(anno_gp_name.attrs[k].astype("str").tolist(), dtype = np.object0)
            k_cate = cate_levels[df[k]]
            df[k] = pd.Categorical(k_cate, categories = cate_levels)
    return df


# recover pandas category attribut

def Py_read_hdf5(file = None):
    if file is None:
        raise OSError('No such file or directory')
    h5 = h5py.File(name=file, mode='r')
    try:
        # normdata
        normcounts = h5_to_matrix(gp_name=h5['normData'])
            # scaledata
        scaledata = h5_to_matrix(gp_name=h5['scaleData'])
            # annotation
        anno = h5['annotation']
        obs_df = h5_to_df(anno_gp_name=anno['observes'])
        sdvar_df = h5_to_df(anno_gp_name=anno['variables/sdVariables'])
        ndvar_df = h5_to_df(anno_gp_name=anno['variables/ndVariables'])
        # create the anndata
        adata = anndata.AnnData(X=scaledata, obs=obs_df, var=sdvar_df)
        adata0 = anndata.AnnData(X=normcounts, obs=obs_df, var=ndvar_df)
        adata.raw = adata0
        # add raw counts
        if 'rawData' in h5.keys():
            counts = h5_to_matrix(gp_name=h5['rawData'])
            adata.layers['counts'] = counts
        # dimension
        dimR = h5['dimReduction']
        for k in dimR.keys():
            X_k = "X_" + k.lower()
            adata.obsm[X_k] = dimR[k][()].astype(np.float32)
        # graphs :
        graphs = h5['graphs']
        neig = {"knn":"distances", "snn":"connectivities"}
        neig_dict ={}
        for g in graphs.keys():
            neig_dict[neig[g]]= h5_to_matrix(graphs[g])
        adata.uns["neighbors"] = neig_dict
    except Exception as e:
        print('Error:',e)
    finally:
        h5.close()
    return adata



# test py read
''' 

'''