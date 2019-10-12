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
import re
import os

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


# metadata to h5
# colors
def meta_to_h5(meta, h5, gp_name=None):
    cr = re.compile('.*('+gp_name+')')
    meta_c = [m.group(0) for l in [l0 for l0 in meta.keys()] for m in [cr.search(l)] if m]
    if len(meta_c)>0:
        h5_gp = h5.create_group('metadata')
        h5_gp_name = h5_gp.create_group(gp_name)
        for c in meta_c:
            h5_gp_name.create_dataset(c, data=meta[c].astype(h5py.special_dtype(vlen=str)))

def adata_to_h5(adata=None, h5=None):
    h5_group={'rawData', 'normData','scaleData','annotation', 'dimReduction', 'graphs', 'metadata'}
    dt = [np.ndarray, scipy.sparse.csr.csr_matrix]
    for k in h5_group:
        if k == 'rawData':
            if 'counts' in adata.uns.keys():
                if type(adata.uns['counts']) in dt:
                    matrix_to_h5(mat=adata.uns['counts'], h5_gp=h5, gp_name=k)
        if k == 'normData':
            ndata = adata.raw.X
            if not ndata is None:
                if type(ndata) in dt:
                    matrix_to_h5(mat=ndata, h5_gp=h5, gp_name=k)
        if k == 'scaleData':
            sdata = adata.X
            if not ndata is None:
                if type(ndata) in dt:
                    matrix_to_h5(mat=sdata, h5_gp=h5, gp_name=k)
        if k == 'annotation':
            annotation = h5.create_group(k)
            if not adata.obs is None:
                df_to_h5(df=adata.obs, h5_anno=annotation, anno_dataset='observes')
            if not adata.var is None:
                df_to_h5(df=adata.var, h5_anno=annotation, anno_dataset='sdVariables')
            if not adata.raw.var is None:
                df_to_h5(df=adata.raw.var, h5_anno=annotation, anno_dataset='ndVariables')
        if k == 'dimReduction':
            if len(adata.obsm.keys()) > 0:
                dimReduction=h5.create_group(k)
                for k in [k for k in adata.obsm.keys()]:
                    K = re.sub("^.*_", "", k).upper()
                    dimReduction.create_dataset(K, data=adata.obsm[k], dtype=np.float32)
        if k == 'graphs':
            if 'neighbors' in adata.uns.keys():
                graphs = h5.create_group(k)
                gra_dict={"distances":"knn", "connectivities":"snn"}
                graph_data = adata.uns['neighbors']
                # -- save the neighbor graphs
                for g in gra_dict.keys():
                    matrix_to_h5(mat=graph_data[g], h5_gp=graphs, gp_name=gra_dict[g])
        if k == 'metadata':
            meta_to_h5(meta=adata.uns, h5=h5, gp_name='colors')
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
        adata_to_h5(adata=adata, h5=h5)
    except Exception as e:
        print('Error:',e)
    finally:
        h5.close()
    return




############ read the hdf5 file

# h5 to the scipy sparse csr format and ndarray
def h5_to_matrix(gp_name):
    if isinstance(gp_name.attrs['datatype'], str):
        if gp_name.attrs['datatype'] == 'SparseMatrix':
            x = gp_name["values"][()].astype(np.float32)
            indices = gp_name["indices"][()]
            indptr = gp_name["indptr"][()]
            shapes = gp_name["dims"][()]
            mat = sparse.csr_matrix((x,indices,indptr),shape=shapes, dtype = np.float32)
        elif gp_name.attrs['datatype'] == 'Array':
            mat = gp_name['matrix'][()].astype(np.float32)
    elif isinstance(gp_name.attrs['datatype'], np.ndarray):
        if gp_name.attrs['datatype'].astype('str').item() =='SparseMatrix':
            x = gp_name["values"][()].astype(np.float32)
            indices = gp_name["indices"][()]
            indptr = gp_name["indptr"][()]
            shapes = gp_name["dims"][()]
            mat = sparse.csr_matrix((x,indices,indptr),shape=shapes, dtype = np.float32)
        elif gp_name.attrs['datatype'].astype('str').item() == 'Array':
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
            k_cate = cate_levels[df[k].astype(np.int64)]
            df[k] = pd.Categorical(k_cate, categories = cate_levels)
    return df


def h5_to_adata(h5 = None):
    if 'normData' in h5.keys():
        normcounts = h5_to_matrix(gp_name=h5['normData'])
        m1 = 1
    if 'scaleData' in h5.keys():
        scaledata = h5_to_matrix(gp_name=h5['scaleData'])
        m2 = 2
    if 'annotation' in h5.keys():
        anno = h5['annotation']
        if 'observes' in anno.keys():
            obs_df = h5_to_df(anno_gp_name=anno['observes'])
            m3 = 3
        if 'sdVariables' in anno.keys():
            sdvar_df = h5_to_df(anno_gp_name=anno['sdVariables'])
            m4 = 4
        if 'ndVariables' in anno.keys():
            ndvar_df = h5_to_df(anno_gp_name=anno['ndVariables'])
            m5 = 5
        if m1+m2+m3+m4+m5 == 15:
            adata = anndata.AnnData(X=scaledata, obs=obs_df, var=sdvar_df)
            adata0 = anndata.AnnData(X=normcounts, obs=obs_df, var=ndvar_df)
            adata.raw = adata0
            if 'rawData' in h5.keys():
                adata.uns['counts'] = h5_to_matrix(gp_name=h5['rawData'])
            if 'dimReduction' in h5.keys():
                dimR = h5['dimReduction']
                for k in dimR.keys():
                    X_k = "X_" + k.lower()
                    adata.obsm[X_k] = dimR[k][()].astype(np.float32)
            if 'graphs' in h5.keys():
                graphs = h5['graphs']
                neig = {"knn":"distances", "snn":"connectivities"}
                neig_dict ={}
                for g in graphs.keys():
                    neig_dict[neig[g]]= h5_to_matrix(graphs[g])
                adata.uns["neighbors"] = neig_dict
            if 'metadata' in h5.keys():
                meta = h5['metadata/colors']
                for k in meta.keys():
                    adata.uns[k] = np.array(meta[k][()].astype("str").tolist(), dtype =np.object0)  
        else:
            raise OSError("Dataset lacks the necessary information")
    return adata
        
# recover pandas category attribut

def Py_read_hdf5(file = None):
    if file is None:
        raise OSError('No such file or directory')
    h5 = h5py.File(name=file, mode='r')
    try:
        adata = h5_to_adata(h5=h5)
    except Exception as e:
        print('Error:',e)
    finally:
        h5.close()
    return adata