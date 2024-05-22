import warnings
warnings.filterwarnings("ignore")
from anndata import AnnData
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import copy

import os
import cv2
import skimage as ski
import matplotlib.pyplot as plt
from PIL import Image
from PIL import ImageFile
ImageFile.LOAD_TRUNCATED_IMAGES = True
Image.MAX_IMAGE_PIXELS = None

import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
import stack_data

import time
from st_processing import dis_two_point,get_chip_coord,get_distance_from_point_to_line,get_HE_mask,get_piexl_coods,get_piexl_coods_Mchip,getLinearEquation_new
from st_processing import rgb2gray,adata_flip,adata_Srotate,trans_pN,make_uniq_list
from st_processing import get_adata_STARsolo,get_adata_STARsolo_10x
from st_processing import load_image,dis_two_point
#bin 14:10um diameter
#bin 50:25um diameter
#bin 140:100um diameter
color_list=['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf',
           '#a6cee3','#fdbf6f','#b2df8a','#fb9a99','#cab2d6','#ffff00','#bd4963','#c5c5fe','#f57682','#6bd76c','#0026bf','#ab4c4c','#116400','#db654e']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.sans-serif'] = 'Arial'




def adata_filter_norm_process(adata_raw,min_cells_num=10,min_genes_num=100,mt_num=20,n_top_genes=12000,use_variable=True,
                              in_tissue=True,image=True,n_clusters=8,scale=False): 
    import time
    start_time = time.time()
    adata=adata_raw.copy()
    if in_tissue:
        adata = adata[adata.obs["in_tissue"] == '1' ]
    if image:
        adata.uns['crop_coord']= [min(adata.obsm['spatial'][:,0])-100, max(adata.obsm['spatial'][:,0])+100,
               min(adata.obsm['spatial'][:,1])-100, max(adata.obsm['spatial'][:,1])+100]

    else:
         adata.uns['crop_coord']= [min(adata.obsm['spatial'][:,0])-3, max(adata.obsm['spatial'][:,0])+3,
               min(adata.obsm['spatial'][:,1])-3, max(adata.obsm['spatial'][:,1])+3]
    if scale:
        sc.pp.filter_cells(adata, min_genes=min_genes_num)
        sc.pp.filter_genes(adata, min_cells=min_cells_num)
        adata = adata[adata.obs["pct_counts_mt"] < mt_num]
        
        adata.layers['raw'] = adata.X.copy()
        sc.pp.normalize_total(adata,target_sum=1e4)
        adata.layers['norm'] = adata.X.copy()
        sc.pp.log1p(adata)
        adata.raw = adata.copy()
        sc.pp.highly_variable_genes(adata,n_top_genes=n_top_genes) 
        adata=adata[:,adata.var.highly_variable]#adata.raw.to_adata()
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
        sc.pp.scale(adata, max_value=10)
        
    else:
        sc.pp.filter_cells(adata, min_genes=min_genes_num)
        sc.pp.filter_genes(adata, min_cells=min_cells_num)
        adata = adata[adata.obs["pct_counts_mt"] < mt_num]
        adata.layers['raw'] = adata.X.copy()
        sc.pp.normalize_total(adata,target_sum=1e4)
        adata.raw = adata.copy()
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata,n_top_genes=n_top_genes) 

    sc.tl.pca(adata, svd_solver="arpack", use_highly_variable=use_variable)
    sc.pp.neighbors(adata,  use_rep='X_pca', n_neighbors=50)
    sc.tl.umap(adata)
    sc.tl.tsne(adata)
    sc.tl.leiden(adata)
    sc.tl.louvain(adata, resolution=1)
    X_pca = adata.obsm['X_pca'] 
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(X_pca) 
    adata.obs['kmeans'] = kmeans.labels_.astype(str)
    
    print('Sample:  '+list(adata.uns['spatial'].keys())[0])
    print(f"#genes after  filter: {adata.n_vars}")
    print(f"#Spot after  filter: {adata.n_obs}")
    print(f"#Median UMIs : {round( np.median(adata.obs['total_counts']),0)}")
    print(f"#Median Genes : {round(np.median(adata.obs['n_genes_by_counts']),0)}")
    print(f"#Mean UMIs : {round( np.mean(adata.obs['total_counts']),0)}")
    print(f"#Mean Genes : {round(np.mean(adata.obs['n_genes_by_counts']),0)}")
    end_time = time.time()
    print("Reading runtime：", round(end_time - start_time,0), "s")
    print('-')
    return adata
def interpolation_M9_150(adata,sample='id',rr=0,cc=0):
    import scipy.sparse as sp
    import random
    import time
    raw_X=adata.X.A.copy()
    adata.obs['Pixel_col']=adata.obsm['spatial'][:,0]
    adata.obs['Pixel_row']=adata.obsm['spatial'][:,1]
    dis=adata.uns['spatial'][sample]['scalefactors']['spot_diameter_fullres']*3
    df_plot=adata[adata.obs["in_tissue"] == '1' ].obs.loc[:,['Spot_col','total_counts']].groupby(['Spot_col']).median()#.plot()
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),' - Start col adjust')
    for c in df_plot.index[np.min(df_plot.index)+1:np.max([i for i in df_plot.index if i <451])-1]:
        if (c in df_plot.index) & (c-1 in df_plot.index) & (c+1 in df_plot.index):
            if df_plot.loc[c,'total_counts'] < (df_plot.loc[c-1,'total_counts']+df_plot.loc[c+1,'total_counts'])/3:
                for i in adata.obs.index[adata.obs['Spot_col']==c]:
                    col=adata.obs.loc[i,'Pixel_col']
                    row=adata.obs.loc[i,'Pixel_row']
                    tmp_X=raw_X[(adata.obs['sample']!='add')&(adata.obs['Pixel_col']<(col+dis)) & (adata.obs['Pixel_col']>(col-dis)) & 
                               (adata.obs['Pixel_row']<(row+dis)) & (adata.obs['Pixel_row']>(row-dis)) & (adata.obs['Spot_col']!=c)]
                    if tmp_X.shape[0] > 0:
                        # pos_tmp=adata.obs.loc[(adata.obs['sample']!='add')&(adata.obs['Pixel_col']<(col+dis)) & (adata.obs['Pixel_col']>(col-dis)) & 
                        #        (adata.obs['Pixel_row']<(row+dis)) & (adata.obs['Pixel_row']>(row-dis)) & (adata.obs['Spot_col']!=c),['Pixel_row','Pixel_col']]
                        # dis_tmp=[dis_two_point((pos_tmp.loc[p,'Pixel_row'],pos_tmp.loc[p,'Pixel_col']),(row,col)) for p in pos_tmp.index]
                        tmp_X=tmp_X[random.randint(0, tmp_X.shape[0]-1),:]#tmp_X[dis_tmp.index(min(dis_tmp)),:]
                        raw_X[adata.obs.index==i]=tmp_X
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),' - Start row adjust')
    df_plot=adata[adata.obs["in_tissue"] == '1' ].obs.loc[:,['Spot_row','total_counts']].groupby(['Spot_row']).median()#.plot()
    for r in df_plot.index[np.min(df_plot.index)+1:np.max([i for i in df_plot.index if i <451])-1]:
        if (r in df_plot.index) & (r-1 in df_plot.index) & (r+1 in df_plot.index):
            if df_plot.loc[r,'total_counts'] < (df_plot.loc[r-1,'total_counts']+df_plot.loc[r+1,'total_counts'])/3:
                for i in adata.obs.index[adata.obs['Spot_row']==r]:
                    col=adata.obs.loc[i,'Pixel_col']
                    row=adata.obs.loc[i,'Pixel_row']
                    tmp_X=raw_X[(adata.obs['sample']!='add')&(adata.obs['Pixel_col']<(col+dis)) & (adata.obs['Pixel_col']>(col-dis)) & 
                               (adata.obs['Pixel_row']<(row+dis)) & (adata.obs['Pixel_row']>(row-dis)) & (adata.obs['Spot_row']!=r)]
                    if tmp_X.shape[0] > 0:
                        # pos_tmp=adata.obs.loc[(adata.obs['sample']!='add')&(adata.obs['Pixel_col']<(col+dis)) & (adata.obs['Pixel_col']>(col-dis)) & 
                        #        (adata.obs['Pixel_row']<(row+dis)) & (adata.obs['Pixel_row']>(row-dis)) & (adata.obs['Spot_row']!=r),['Pixel_row','Pixel_col']]
                        # dis_tmp=[dis_two_point((pos_tmp.loc[p,'Pixel_row'],pos_tmp.loc[p,'Pixel_col']),(row,col)) for p in pos_tmp.index]
                        tmp_X=tmp_X[random.randint(0, tmp_X.shape[0]-1),:]#tmp_X[dis_tmp.index(min(dis_tmp)),:]
                        raw_X[adata.obs.index==i]=tmp_X
    if (cc>0) & (rr>0):
        for i in adata.obs.index[(adata.obs['Spot_col']==cc) & (adata.obs['Spot_row']<rr) ]:
            col=adata.obs.loc[i,'Pixel_col']
            row=adata.obs.loc[i,'Pixel_row']
            tmp_X=raw_X[(adata.obs['sample']!='add')&(adata.obs['Pixel_col']<(col+dis)) & (adata.obs['Pixel_col']>(col-dis)) & 
                       (adata.obs['Pixel_row']<(row+dis)) & (adata.obs['Pixel_row']>(row-dis)) & (adata.obs['Spot_col']!=cc)]
            if tmp_X.shape[0] > 0:
                # pos_tmp=adata.obs.loc[(adata.obs['sample']!='add')&(adata.obs['Pixel_col']<(col+dis)) & (adata.obs['Pixel_col']>(col-dis)) & 
                #        (adata.obs['Pixel_row']<(row+dis)) & (adata.obs['Pixel_row']>(row-dis)) & (adata.obs['Spot_col']!=cc),['Pixel_row','Pixel_col']]
                # dis_tmp=[dis_two_point((pos_tmp.loc[p,'Pixel_row'],pos_tmp.loc[p,'Pixel_col']),(row,col)) for p in pos_tmp.index]
                tmp_X=tmp_X[random.randint(0, tmp_X.shape[0]-1),:]#tmp_X[dis_tmp.index(min(dis_tmp)),:]
                raw_X[adata.obs.index==i]=tmp_X
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),' - Start col interpolation')
    dis=adata.uns['spatial'][sample]['scalefactors']['spot_diameter_fullres']*7.5
    l=list(adata.obs['Spot_col'])
    adata.obs['sample_2']=adata.obs['sample']
    for c in [451,457,452,456,453,455,454,458,464,459,463,460,462,461]:
        adata.obs['sample_2'][adata.obs['Spot_col']==c]=='add2'
        if (c in l):
            for i in adata.obs.index[adata.obs['Spot_col']==c]:
                col=adata.obs.loc[i,'Pixel_col']
                row=adata.obs.loc[i,'Pixel_row']
                if adata.obs.loc[i,'Spot_row']>450:
                    tmp_X=raw_X[(adata.obs['sample_2']!='add')&(adata.obs['Pixel_col']<(col+dis*2)) & (adata.obs['Pixel_col']>(col-dis*2)) & 
                               (adata.obs['Pixel_row']<(row+dis*2)) & (adata.obs['Pixel_row']>(row-dis*2)) & (adata.obs['Spot_col']!=c)]
                    # pos_tmp=adata.obs.loc[(adata.obs['sample_2']!='add')&(adata.obs['Pixel_col']<(col+dis*2)) & (adata.obs['Pixel_col']>(col-dis*2)) & 
                    #            (adata.obs['Pixel_row']<(row+dis*2)) & (adata.obs['Pixel_row']>(row-dis*2)) & (adata.obs['Spot_col']!=c),['Pixel_row','Pixel_col']]
                else:
                    tmp_X=raw_X[(adata.obs['sample_2']!='add')&(adata.obs['Pixel_col']<(col+dis)) & (adata.obs['Pixel_col']>(col-dis)) & 
                               (adata.obs['Pixel_row']<(row+dis)) & (adata.obs['Pixel_row']>(row-dis)) & (adata.obs['Spot_col']!=c)]
                    # pos_tmp=adata.obs.loc[(adata.obs['sample_2']!='add')&(adata.obs['Pixel_col']<(col+dis)) & (adata.obs['Pixel_col']>(col-dis)) & 
                    #            (adata.obs['Pixel_row']<(row+dis)) & (adata.obs['Pixel_row']>(row-dis)) & (adata.obs['Spot_col']!=c),['Pixel_row','Pixel_col']]
                if tmp_X.shape[0] > 0:
                    #dis_tmp=[dis_two_point((pos_tmp.loc[p,'Pixel_row'],pos_tmp.loc[p,'Pixel_col']),(row,col)) for p in pos_tmp.index]
                    tmp_X=tmp_X[random.randint(0, tmp_X.shape[0]-1),:]#tmp_X[dis_tmp.index(min(dis_tmp)),:]
                    raw_X[adata.obs.index==i]=tmp_X
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),' - Start row interpolation')
    l=list(adata.obs['Spot_row'])
    for r in [451,457,452,456,453,455,454,458,464,459,463,460,462,461]:
        adata.obs['sample_2'][adata.obs['Spot_row']==r]=='add2'
        if (r in l):
            for i in adata.obs.index[adata.obs['Spot_row']==r]:
                col=adata.obs.loc[i,'Pixel_col']
                row=adata.obs.loc[i,'Pixel_row']
                if adata.obs.loc[i,'Spot_col']>450:
                    tmp_X=raw_X[(adata.obs['sample_2']!='add')&(adata.obs['Pixel_col']<(col+dis*2)) & (adata.obs['Pixel_col']>(col-dis*2)) & 
                               (adata.obs['Pixel_row']<(row+dis*2)) & (adata.obs['Pixel_row']>(row-dis*2)) & (adata.obs['Spot_row']!=r)]
                    # pos_tmp=adata.obs.loc[(adata.obs['sample_2']!='add')&(adata.obs['Pixel_col']<(col+dis*2)) & (adata.obs['Pixel_col']>(col-dis*2)) & 
                    #            (adata.obs['Pixel_row']<(row+dis*2)) & (adata.obs['Pixel_row']>(row-dis*2)) & (adata.obs['Spot_row']!=r),['Pixel_row','Pixel_col']]
                else:
                    tmp_X=raw_X[(adata.obs['sample_2']!='add')&(adata.obs['Pixel_col']<(col+dis)) & (adata.obs['Pixel_col']>(col-dis)) & 
                               (adata.obs['Pixel_row']<(row+dis)) & (adata.obs['Pixel_row']>(row-dis)) & (adata.obs['Spot_row']!=r)]
                    # pos_tmp=adata.obs.loc[(adata.obs['sample_2']!='add')&(adata.obs['Pixel_col']<(col+dis)) & (adata.obs['Pixel_col']>(col-dis)) & 
                    #            (adata.obs['Pixel_row']<(row+dis)) & (adata.obs['Pixel_row']>(row-dis)) & (adata.obs['Spot_row']!=r),['Pixel_row','Pixel_col']]
                if tmp_X.shape[0] > 0:
                    #dis_tmp=[dis_two_point((pos_tmp.loc[p,'Pixel_row'],pos_tmp.loc[p,'Pixel_col']),(row,col)) for p in pos_tmp.index]
                    tmp_X=tmp_X[random.randint(0, tmp_X.shape[0]-1),:]#tmp_X[dis_tmp.index(min(dis_tmp)),:]
                    raw_X[adata.obs.index==i]=tmp_X
                    
    adata.X=sp.csr_matrix(raw_X)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt","rp","hb"], inplace=True)
    return adata


geneinfo_path='/Public_data/reference/mm10/'


sample='E17-1'
file_path='/MAGIC_seq/Mouse_Embryo_M9_150_20um/data/result_STARsolo/'+sample+'/STARsolo/'
image_file_path='/MAGIC_seq/Mouse_Embryo_M9_150_20um/data/Image/'
Barcode_file_path='/MAGIC_seq/Mouse_Embryo_M9_150_20um/data/Barcode-M9-150-E17.5/'

res_um=20 #点直径
channels_num=450 #通道数
barcode_num=150 #barcode数
reg='reg1' #区域
chip_type='M9'

print(f'sample: {sample}')

HE_point1=(98-5,179-15)
Spot_point1=(161,242)
HE_point2=(4704-12,4933-36)
Spot_point2=(4815,4738)#

line_point_r1c1=(311,377) 
line_point_r1c70=(398,4645)
line_point_r70c1=(4574,322)



adata=get_adata_STARsolo(sample,chip_type,reg,channels_num,barcode_num,res_um=res_um,add_M9_20=True,
                        img_file=True,EM=True,Velocyto=True,soloFeatures='GeneFull',raw=True,species='mouse',
                        file_path=file_path,image_file_path=image_file_path,Barcode_file_path=Barcode_file_path,geneinfo_path=geneinfo_path,
                        HE_point1=HE_point1,Spot_point1=Spot_point1,HE_point2=HE_point2,Spot_point2=Spot_point2,
                        line_point_r1c1=line_point_r1c1,line_point_r1c70=line_point_r1c70,line_point_r70c1=line_point_r70c1)
adata.write('/MAGIC_seq/file/Mouse_Embryo_M9_150_20um/st_adata/adata_'+sample+'.h5ad')
del adata



sample='P0-1'
file_path='/MAGIC_seq/Mouse_Embryo_M9_150_20um/data/result_STARsolo/'+sample+'/STARsolo/'
image_file_path='/MAGIC_seq/Mouse_Embryo_M9_150_20um/data/Image/'
Barcode_file_path='/MAGIC_seq/Mouse_Embryo_M9_150_20um/data/Barcode-M9-150-P04/'

res_um=20 #点直径
channels_num=450 #通道数
barcode_num=150 #barcode数
reg='reg1' #区域
chip_type='M9'

print(f'sample: {sample}')

HE_point1=(124,249)
Spot_point1=(189,334)
HE_point2=(4875,4896)
Spot_point2=(4841,4724)#

line_point_r1c1=(375,424) 
line_point_r1c70=(404,4693)
line_point_r70c1=(4631,375)
adata=get_adata_STARsolo(sample,chip_type,reg,channels_num,barcode_num,res_um=res_um,add_M9_20=True,
                        img_file=True,EM=True,Velocyto=True,soloFeatures='GeneFull',raw=True,species='mouse',
                        file_path=file_path,image_file_path=image_file_path,Barcode_file_path=Barcode_file_path,geneinfo_path=geneinfo_path,
                        HE_point1=HE_point1,Spot_point1=Spot_point1,HE_point2=HE_point2,Spot_point2=Spot_point2,
                        line_point_r1c1=line_point_r1c1,line_point_r1c70=line_point_r1c70,line_point_r70c1=line_point_r70c1)
adata.write('/MAGIC_seq/file/Mouse_Embryo_M9_150_20um/st_adata/adata_'+sample+'.h5ad')
del adata

sample='P4-1'
file_path='/MAGIC_seq/Mouse_Embryo_M9_150_20um/data/result_STARsolo/'+sample+'/STARsolo/'
image_file_path='/MAGIC_seq/Mouse_Embryo_M9_150_20um/data/Image/'
Barcode_file_path='/MAGIC_seq/Mouse_Embryo_M9_150_20um/data/Barcode-M9-150-P04/'

res_um=20 #点直径
channels_num=450 #通道数
barcode_num=150 #barcode数
reg='reg1' #区域
chip_type='M9'

print(f'sample: {sample}')

HE_point1=(162,195)
Spot_point1=(115,371)
HE_point2=(4780,4861)
Spot_point2=(4461,4682)#

line_point_r1c1=(254,501) 
line_point_r1c70=(289,4527)
line_point_r70c1=(4283,519)

adata=get_adata_STARsolo(sample,chip_type,reg,channels_num,barcode_num,res_um=res_um,add_M9_20=True,
                        img_file=True,EM=True,Velocyto=True,soloFeatures='GeneFull',raw=True,species='mouse',
                        file_path=file_path,image_file_path=image_file_path,Barcode_file_path=Barcode_file_path,geneinfo_path=geneinfo_path,
                        HE_point1=HE_point1,Spot_point1=Spot_point1,HE_point2=HE_point2,Spot_point2=Spot_point2,
                        line_point_r1c1=line_point_r1c1,line_point_r1c70=line_point_r1c70,line_point_r70c1=line_point_r70c1)
adata.write('/MAGIC_seq/file/Mouse_Embryo_M9_150_20um/st_adata/adata_'+sample+'.h5ad')
del adata

sample='E17-1'
adata_E17=sc.read('/MAGIC_seq/file/Mouse_Embryo_M9_150_20um/st_adata/adata_'+sample+'.h5ad')
adata_E17.var['gene_id']=adata_E17.var_names
adata_E17.var_names=make_uniq_list(adata_E17.var['Symbol'])
print(sample)
print('Median UMIs:',round(np.median(adata_E17.obs[(adata_E17.obs['sample']!='add') & (adata_E17.obs['in_tissue']!='0')]['total_counts']),0))
print('Median Genes:',round(np.median(adata_E17.obs[(adata_E17.obs['sample']!='add') & (adata_E17.obs['in_tissue']!='0')]['n_genes_by_counts']),0))
print('In tissue(%):',round(100-100*np.sum(adata_E17.obs['total_counts'][(adata_E17.obs['sample']!='add') & (adata_E17.obs['in_tissue']=='0')])/np.sum(adata_E17.obs['total_counts']),0))
adata_E17 = adata_E17[adata_E17.obs["in_tissue"] == '1' ]
sc.pp.filter_genes(adata_E17, min_cells=20)
adata_E17=interpolation_M9_150(adata_E17,sample=sample)
adata_E17.write('/MAGIC_seq/file/Mouse_Embryo_M9_150_20um/st_adata/adata_'+sample+'_add.h5ad')

adata_E17_filtered=adata_filter_norm_process(adata_E17,min_cells_num=20,min_genes_num=100,mt_num=30,n_top_genes=16000,use_variable=True,
                                     in_tissue=True,image=True,n_clusters=16)
adata_E17_filtered.write('/MAGIC_seq/file/Mouse_Embryo_M9_150_20um/st_adata/adata_'+sample+'_filtered.h5ad')


sample='P0-1'
adata_P0=sc.read('/MAGIC_seq/file/Mouse_Embryo_M9_150_20um/st_adata/adata_'+sample+'.h5ad')
adata_P0.var['gene_id']=adata_P0.var_names
adata_P0.var_names=make_uniq_list(adata_P0.var['Symbol'])
print(sample)
print('Median UMIs:',round(np.median(adata_P0.obs[(adata_P0.obs['sample']!='add') & (adata_P0.obs['in_tissue']!='0')]['total_counts']),0))
print('Median Genes:',round(np.median(adata_P0.obs[(adata_P0.obs['sample']!='add') & (adata_P0.obs['in_tissue']!='0')]['n_genes_by_counts']),0))
print('In tissue(%):',round(100-100*np.sum(adata_P0.obs['total_counts'][(adata_P0.obs['sample']!='add') & (adata_P0.obs['in_tissue']=='0')])/np.sum(adata_P0.obs['total_counts']),0))
adata_P0 = adata_P0[adata_P0.obs["in_tissue"] == '1' ]
sc.pp.filter_genes(adata_P0, min_cells=20)
adata_P0=interpolation_M9_150(adata_P0,sample=sample)
adata_P0.write('/MAGIC_seq/file/Mouse_Embryo_M9_150_20um/st_adata/adata_'+sample+'_add.h5ad')
adata_P0_filtered=adata_filter_norm_process(adata_P0,min_cells_num=20,min_genes_num=100,mt_num=30,n_top_genes=16000,use_variable=True,
                                     in_tissue=True,image=True,n_clusters=16)
adata_P0_filtered.write('/MAGIC_seq/file/Mouse_Embryo_M9_150_20um/st_adata/adata_'+sample+'_filtered.h5ad')


sample='P4-1'
adata_P4=sc.read('/MAGIC_seq/file/Mouse_Embryo_M9_150_20um/st_adata/adata_'+sample+'.h5ad')
adata_P4.var['gene_id']=adata_P4.var_names
adata_P4.var_names=make_uniq_list(adata_P4.var['Symbol'])
print(sample)
print('Median UMIs:',round(np.median(adata_P4.obs[(adata_P4.obs['sample']!='add') & (adata_P4.obs['in_tissue']!='0')]['total_counts']),0))
print('Median Genes:',round(np.median(adata_P4.obs[(adata_P4.obs['sample']!='add') & (adata_P4.obs['in_tissue']!='0')]['n_genes_by_counts']),0))
print('In tissue(%):',round(100-100*np.sum(adata_P4.obs['total_counts'][(adata_P4.obs['sample']!='add') & (adata_P4.obs['in_tissue']=='0')])/np.sum(adata_P4.obs['total_counts']),0))
adata_P4 = adata_P4[adata_P4.obs["in_tissue"] == '1' ]
sc.pp.filter_genes(adata_P4, min_cells=20)
adata_P4=interpolation_M9_150(adata_P4,sample=sample,cc=301,rr=80)
adata_P4.write('/MAGIC_seq/file/Mouse_Embryo_M9_150_20um/st_adata/adata_'+sample+'_add.h5ad')

adata_P4_filtered=adata_filter_norm_process(adata_P4,min_cells_num=20,min_genes_num=100,mt_num=30,n_top_genes=16000,use_variable=True,
                                     in_tissue=True,image=True,n_clusters=16)
adata_P4_filtered.write('/MAGIC_seq/file/Mouse_Embryo_M9_150_20um/st_adata/adata_'+sample+'_filtered.h5ad')
