def read_10X_visium(file_path,count_file):
    import scanpy as sc
    import copy
    adata=sc.read_visium(file_path,count_file=count_file)
    adata.obsm['spatial']=adata.obsm['spatial'].astype(int)
    adata.obs['in_tissue']=adata.obs['in_tissue'].astype(str)
    sample=list(adata.uns['spatial'].keys())[0]
    adata.var_names_make_unique()
    adata.layers['raw']=copy.deepcopy(adata.X)
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    adata.var["hb"] = adata.var_names.str.startswith("Hb-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt","hb"], inplace=True)
    print(f'{sample}:')
    print('Median total_counts',np.median(adata[adata.obs["in_tissue"] == '1'  ].obs['total_counts']))
    print('Median n_genes_by_counts',np.median(adata[adata.obs["in_tissue"] == '1'  ].obs['n_genes_by_counts']))
    print('-')
    return adata
