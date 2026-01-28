# useful functions for single cell analysis

def BarplotOrig(which_var, adata, var='clusters', height=3, color = False):
    import pandas as pd
    import scanpy as sc
    plotdata = pd.crosstab(adata.obs[var], adata.obs[which_var], normalize='index') * 100
    if 'category' in plotdata.index.dtype.name:
        plotdata.index.reorder_categories(adata.obs[var].cat.categories[::-1])

    if not color:
        ax1 = plotdata.plot.barh(stacked = True, edgecolor = 'none', zorder = 3, figsize = (6,height),
                                 fontsize = 14, grid = False)
    else:
        ax1 = plotdata.plot.barh(stacked = True, edgecolor = 'none', zorder = 3, figsize = (6,height),
                                 fontsize = 14, grid = False, color = color)
    ax1.set_title(which_var+' %')
    ax1.set_ylabel(var)
    horiz_offset = 1
    vert_offset = 1
    ax1 = ax1.legend(bbox_to_anchor = (horiz_offset, vert_offset))
    ax1.figure.savefig(str(sc.settings.figdir)+'/barplot_'+var+'_proportions_'+which_var+'.pdf', bbox_inches='tight',
                       dpi=300, orientation='landscape', format= 'pdf')
    # ax1.figure.savefig(str(sc.settings.figdir)+'/barplot_'+var+'_proportions_'+which_var+'.pdf', bbox_inches='tight',
    #                    dpi=300, orientation='landscape', format= 'pdf', optimize=True)
    

def Barplot(which_var, adata, var='clusters', figsize = (6,3), color = None, fontsize = 14, grid = None, var_order = None, group_order = None, **kwargs):
    import pandas as pd
    import scanpy as sc
    plotdata = pd.crosstab(adata.obs[var], adata.obs[which_var], normalize='index') * 100
    if var_order: 
        plotdata = plotdata.loc[var_order]
    if group_order:
        plotdata = plotdata[group_order]
    if kwargs.get('sort_by'):
        if not kwargs.get('ascending'):
            kwargs['ascending'] = False
        plotdata = plotdata.sort_values(by = kwargs.get('sort_by'), ascending=kwargs.get('ascending'))
    if 'category' in plotdata.index.dtype.name:
        plotdata.index.reorder_categories(adata.obs[var].cat.categories[::-1])
    if not color:
        ax1 = plotdata.plot.barh(stacked = True, edgecolor = 'none', zorder = 3, figsize = figsize,
                                 fontsize = fontsize, grid = grid);
    else:
        ax1 = plotdata.plot.barh(stacked = True, edgecolor = 'none', zorder = 3, figsize = figsize,
                                 fontsize = fontsize, grid = grid, color = color);
    ax1.set_title(which_var+' %');
    ax1.set_ylabel(var);
    return ax1

#     ax1.figure.savefig(str(sc.settings.figdir)+'/barplot_'+var+'_proportions_'+which_var+'.pdf', bbox_inches='tight',
#                        dpi=300, orientation='landscape', format= 'pdf')

def dge2df(obj, group, QVAL_THRESH, LOGFC_THRESH):
    import pandas as pd
    # group = 'thy2'
    # QVAL_THRESH = 0.2
    # LOGFC_THRESH = 1
    t1 = pd.DataFrame(obj.uns['rank_genes_groups']['names'])[group].rename('names')
    t2 = pd.DataFrame(obj.uns['rank_genes_groups']['scores'])[group].rename('scores')
    t3 = pd.DataFrame(obj.uns['rank_genes_groups']['pvals'])[group].rename('pvals')
    t4 = pd.DataFrame(obj.uns['rank_genes_groups']['pvals_adj'])[group].rename('pvals_adj')
    t5 = pd.DataFrame(obj.uns['rank_genes_groups']['logfoldchanges'])[group].rename('logfoldchanges')
    res = pd.concat([t1,t2,t3,t4,t5], axis=1)
    outputFile = res[(res['pvals_adj'] < QVAL_THRESH) & (res['logfoldchanges'] > LOGFC_THRESH)].sort_values('scores', ascending=False)
    # print(f'outputFile size: {outputFile.shape}')
    return outputFile

def makeOutoutFolders():
    import pathlib as pt
    # make output path
    outFigurs = './figures/'
    outResults = './results/'
    pt.Path(outFigurs).mkdir(parents=True, exist_ok=True)
    pt.Path(outResults).mkdir(parents=True, exist_ok=True)
    return outFigurs, outResults

def nbInit(outFigurs):
    import scanpy as sc
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import session_info
    # scanpy settings
    sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.settings.figdir = outFigurs
    sc.settings.set_figure_params(dpi=80)  # low dpi (dots per inch) yields small inline figures
    mpl.rcParams['pdf.fonttype'] = 42 #to export full text strings to illustrator
    mpl.rcParams['ps.fonttype'] = 42
    font = {'size'   : 10}
    mpl.rc('font', **font)
    mpl.rc('xtick', labelsize=10) 
    mpl.rc('ytick', labelsize=10) 
    axes = {'labelsize': 10,
            'titlesize': 12}
    mpl.rc('axes', **axes)
    np.random.seed(0)
    
def summaryTable(adata, ll, return_df=False, **save_kwargs):
    import pandas as pd
    a = set(ll)
    if not ll: #add this columns if a is empty
        b = set(['donor','pcw','sex','gender']) 
    else: b=[]
    tempAnn = adata.obs[list(set(adata.obs_keys()).intersection(list(a.union(b))))]
    if tempAnn.columns.isin(['pcw']).any():
        try: tempAnn['pcw'] = tempAnn['pcw'].values.astype(int)
        except ValueError: pass
        sampleTable = pd.DataFrame(tempAnn.value_counts().sort_index(level='pcw', ascending=True).rename('cells').reset_index())
    else:
        sampleTable = pd.DataFrame(tempAnn.value_counts().rename('cells').reset_index())
    print(f"[INFO] - total cells: {sampleTable['cells'].sum()}")
    if save_kwargs:
        sampleTable.to_csv(save_kwargs.pop('save'), **save_kwargs)
    if return_df:
        return sampleTable
    else:
        return(sampleTable.style)

def getExpPercGroup(adata, groups, min_cells = 20,  min_pct = 0.1, show_aggs=True):
    '''return expressed cells percentage and number per group, output genes and groups
    sf.getExpPercGroup(adata[adata.obs['celltype'].isin(['thy_TH_processing', 'thy_Lumen-forming'])], ['age_group','celltype','karyotype'])'''
    import numpy as np
    import pandas as pd
    df = adata.obs[groups].agg('_'.join, axis=1).rename('agg_col').reset_index()
    agg_dict = df.groupby('agg_col')['index'].apply(list).to_dict()
    perc_dict = {}
    exp_dict = {}
    for k,v in agg_dict.items():
        siz = np.size(adata[v,:].X, axis=0)
        exp = np.sum(adata[v,:].X > 0, axis=0).A1 
        perc_dict[k] = exp/siz
        exp_dict[k] = exp
        
    p = pd.DataFrame.from_dict(perc_dict).set_index(adata.var_names)
    e = pd.DataFrame.from_dict(exp_dict).set_index(adata.var_names)
    HI_GENES = ((e >= min_cells) * (p >= min_pct)).any(axis=1)
    HI_GENES = HI_GENES[HI_GENES].index.tolist()
    print(f'[INFO] {len(agg_dict.keys())} aggs using {groups} at min_cells={min_cells},  min_pct={min_pct}. {len(HI_GENES)} genes')
    if show_aggs:
        display(list(agg_dict.keys()))
    return HI_GENES, list(agg_dict.keys())

def getProteinCodingGenes(adata, group='protein_coding'):
    '''return protein-coding genes in adata
    example: getProteinCodingGenes(adata)'''
    import pandas as pd
    coding_genes = pd.read_csv('/nfs/team292/hm11/endo_GLND/starsolo/data/human_2020A/gene_symbol_type.tsv', sep='\t', header=None)
    cgenes = coding_genes.loc[coding_genes[2] == 'protein_coding',1].tolist()
    cgenesInData_logic = adata.var_names.isin(cgenes)
    print(f'[INFO] coding genes found in adata {cgenesInData_logic.sum()} ({len(cgenes)} total protein-coding genes)')
    return adata.var.loc[cgenesInData_logic,:].index.tolist()

def makeColors(cmap='flare', N=5, showCols=None, **kwargs):
    from matplotlib import pyplot as plt
    import matplotlib as mpl
    import seaborn as sns
    pallcol = sns.color_palette(cmap, N, **kwargs);
    if showCols:
        sns.palplot(pallcol);
    return list(map(mpl.colors.rgb2hex, pallcol)), pallcol

# readExcelColors
def readExcelColors(colname: str, return_df=False):
    '''return colors from excel file in Cat order by adata'''
    import pandas as pd
    colors = pd.read_excel('/nfs/team292/hm11/endo_GLND/starsolo/GLND/Thyroid/thyroid_color_code.xlsx')
    tmp = colors.loc[:,[colname,f'{colname}_hex_code']].drop_duplicates().set_index(colname)
    if return_df: return tmp
    else: return tmp.to_dict()[f'{colname}_hex_code']


# # update colors based on DICT for obsColumn
# def updateColors(adata, obsColumn, DICT):
#     '''adata.uns[f'{obsColumn}_colors']'''
#     import pandas as pd
#     inputTbl = pd.DataFrame.from_dict(DICT, orient='index', columns=['color'])
#     if f'{obsColumn}_colors' in adata.uns:
#         tmp = pd.DataFrame({'color'  : adata.uns[f'{obsColumn}_colors']}, index =adata.obs[obsColumn].astype('category').cat.categories)
#         tmp.update(inputTbl)
#         return tmp['color'].str.upper().values.tolist()
#     else:
#         print(f'{obsColumn} not found in adata.uns. Generating adata.uns...')
#         adata.uns[f'{obsColumn}_colors'] = list(adata.obs[obsColumn].astype('category').cat.categories)
#         tmp = pd.DataFrame({'color'  : adata.uns[f'{obsColumn}_colors']}, index =adata.obs[obsColumn].astype('category').cat.categories)
#         tmp.update(inputTbl)
#         return tmp['color'].str.upper().values.tolist()

# update colors based on DICT for obsColumn
def updateColors(adata, obsColumn, DICT):
    '''adata.uns[f'{obsColumn}_colors']'''
    import pandas as pd
    inputTbl = pd.DataFrame.from_dict(DICT, orient='index', columns=['colors'])
    if f'{obsColumn}_colors' in adata.uns:
        # tmp = pd.DataFrame({'color'  : adata.uns[f'{obsColumn}_colors']}, index =adata.obs[obsColumn].astype('category').cat.categories)
        tmp = getCatColors(adata, obsColumn, return_df=True)
        intersect = list(set(inputTbl.index) - set(tmp.index))
        if intersect:
            print(f"[INFO] adding index {intersect} to color list...")
            tmp = pd.concat([tmp,inputTbl.loc[intersect]], axis=0)
        tmp.update(inputTbl)
        return tmp['colors'].str.upper().values.tolist()
    else:
        print(f'{obsColumn} not found in adata.uns. Generating adata.uns...')
        adata.uns[f'{obsColumn}_colors'] = list(adata.obs[obsColumn].astype('category').cat.categories)
#         tmp = pd.DataFrame({'color'  : adata.uns[f'{obsColumn}_colors']}, index =adata.obs[obsColumn].astype('category').cat.categories)
        tmp = getCatColors(adata, obsColumn, return_df=True)
        intersect = list(set(inputTbl.index) - set(tmp.index))
        if intersect:
            print(f"[INFO] adding index {intersect} to color list...")
            tmp = pd.concat([tmp,inputTbl.loc[intersect]], axis=0)
        tmp.update(inputTbl)
        return tmp['colors'].str.upper().values.tolist()

def getCatColors(adata, obsCat, return_df=False):
    '''get colors of categorial column obsCat'''
    import pandas as pd
    if not adata.obs[obsCat].dtype == 'category':
        adata.obs[obsCat] = adata.obs[obsCat].astype('category')
        # print(f'setting {obsCat} as categorial')
    cats = adata.obs[obsCat].cat.categories.tolist()
    cols = adata.uns[f'{obsCat}_colors']
    tmp = pd.DataFrame(list(zip(cats, cols)), columns=[obsCat, 'colors'])
    if return_df:
        return tmp.set_index(obsCat)
    if not return_df:
        return cols

# factor color
def factorColor(adata, obsColumn, cluster, FACTOR):
    from PIL import ImageColor
    import matplotlib
    import pandas as pd 
    import numpy as np
    
    # factor color function
    def updateColor(hex_color, FACTOR):
        rgb = np.array(ImageColor.getrgb(hex_color))
        factor_rgb = (((255 - rgb) * FACTOR) + rgb)/255
        return matplotlib.colors.to_hex(factor_rgb, keep_alpha=True)
    
    # update color
    tmp = pd.DataFrame({obsColumn : adata.obs[obsColumn].cat.categories,
                        'color'  : adata.uns[f'{obsColumn}_colors']})
    cluster = [cluster] if isinstance(cluster, str) else cluster
    hex_new = tmp.loc[tmp[obsColumn].isin(cluster),'color'].apply(lambda c: updateColor(c, FACTOR))
    tmp['color'].update(hex_new)
    return tmp['color'].values

def removeGenes(adata, min_counts=3):
    import pandas as pd
    import scanpy as sc

    # get adata shape 
    shape0 = adata.shape[1]

    # remove stress genes
    df = pd.read_csv('/nfs/team292/hm11/endo_GLND/starsolo/data/stressgenes_vandenBrink2017.csv') #PMID: 28960196
    stress_genes = df['gene'].tolist()
    stress_genes = [ i.upper() for i in stress_genes]
    not_sgenes = [i not in stress_genes for i in adata.var_names]
    adata = adata[:,not_sgenes]
    shape1 = adata.shape[1]

    # remove mito genes
    non_mito_genes = [name for name in adata.var_names if not name.startswith('MT-')]
    adata = adata[:, non_mito_genes]
    shape2 = adata.shape[1]

    # remove ribo genes
    non_ribo_genes = [name for name in adata.var_names if not name.startswith('RP')]
    adata = adata[:, non_ribo_genes]
    shape3 = adata.shape[1]

    # remove lowly expressed genes
    sc.pp.filter_genes(adata, min_counts=min_counts)
    shape4 = adata.shape[1]

    # print msg 
    print(f'removing {shape0 - shape1} stress genes')
    print(f'removing {shape1 - shape2} mito genes')
    print(f'removing {shape2 - shape3} ribo genes')
    print(f'removing {shape3 - shape4} lowly expressed genes [min_counts {min_counts}]')
    print(f'adata shape {adata.shape}')
    return adata

def addChrInfo(df, geneColumn=None, mark=None):
    # /add chr information. See path in the functions/#
    #features.tsv file from the cellranger output, and then add the chromosome information from the gtf file
    # dependent packages
    import pandas as pd
    # load chr info
    path = '/nfs/team292/hm11/endo_GLND/starsolo/data/fThy_geneAnnotation_GRCh38-2020-A.csv'
    chrinfo = pd.read_csv(path)
    chrinfo = chrinfo.dropna().drop_duplicates(subset='geneSym').rename(columns={'geneSym':'gene'})[['gene','chr','ensID']]
    chrinfo = chrinfo.astype(str)
    # add mark
    if mark != None:
        chrinfo.loc[chrinfo['chr'] != mark, 'chr'] = 'other'
        chrinfo.loc[chrinfo['chr'] == mark, 'chr'] = mark
    # looks for gene column
    if geneColumn == None: 
        df = df.reset_index(names=['gene'])
    # concat chr information
    df = df.drop(columns=['chr','strand'], errors='ignore').merge(chrinfo, how='left', left_on='gene', right_on='gene')
    # reset mat as the forme omported to the function
    if geneColumn == None:
        df = df.set_index('gene')
    return df

def adjust_box_widths(g, fac):
    """
    Adjust the widths of a seaborn-generated boxplot.
    """
    
    # dependent packages
    import numpy as np
    from matplotlib.patches import PathPatch

    # iterating through Axes instances
    for ax in g.axes:

        # iterating through axes artists:
        for c in ax.get_children():

            # searching for PathPatches
            if isinstance(c, PathPatch):
                # getting current width of box:
                p = c.get_path()
                verts = p.vertices
                verts_sub = verts[:-1]
                xmin = np.min(verts_sub[:, 0])
                xmax = np.max(verts_sub[:, 0])
                xmid = 0.5*(xmin+xmax)
                xhalf = 0.5*(xmax - xmin)

                # setting new width of box
                xmin_new = xmid-fac*xhalf
                xmax_new = xmid+fac*xhalf
                verts_sub[verts_sub[:, 0] == xmin, 0] = xmin_new
                verts_sub[verts_sub[:, 0] == xmax, 0] = xmax_new

                # setting new width of median line
                for l in ax.lines:
                    if np.all(l.get_xdata() == [xmin, xmax]):
                        l.set_xdata([xmin_new, xmax_new])

# filter functions
def filterFM(df, celltype, FDR_THRESH=0.01, LOGFC_THRESH=0.2, sortby='LOGFC', print_bins=None):
    '''filter df and sort by logFC'''
    import pandas as pd
    bins=[0,1.00e-04,1.00e-03,1.00e-02,5.00e-02,1]
    labels=['****','***','**','*','ns']
    df['p_val_adj_sym'] = pd.cut(df['p_val_adj'], bins=bins, right=True, labels=labels)
    if type(celltype) is not list: celltype = [celltype]
    up = df[(df['p_val_adj'] < FDR_THRESH) & (df['avg_log2FC'] > LOGFC_THRESH)  & (df['celltype'].isin(celltype))].sort_values(['celltype','avg_log2FC'],ascending=False)
    dn = df[(df['p_val_adj'] < FDR_THRESH) & (df['avg_log2FC'] < -LOGFC_THRESH) & (df['celltype'].isin(celltype))].sort_values(['celltype','avg_log2FC'],ascending=True)
    if sortby!='LOGFC':
        up = up.sort_values(['celltype','p_val_adj'],ascending=True)
        dn = dn.sort_values(['celltype','p_val_adj'],ascending=True)
    print(f'filtering at FDR_THRESH={FDR_THRESH}, LOGFC_THRESH={LOGFC_THRESH}, length={len(up)}/{len(dn)}/{len(up)+len(dn)}')
    if print_bins:
        print(dict(zip(bins,labels)))
    return up, dn

def filterEdgeR(df, celltype, FDR_THRESH=0.01, LOGFC_THRESH=0.2, sortby='LOGFC', print_bins=None):
    '''filter df and sort by logFC'''
    import pandas as pd
    bins=[0,1.00e-04,1.00e-03,1.00e-02,5.00e-02,1]
    labels=['****','***','**','*','ns']
    df['FDR_sym'] = pd.cut(df['FDR'], bins=bins, right=True, labels=labels)
    if type(celltype) is not list: celltype = [celltype]
    up = df[(df['FDR'] < FDR_THRESH) & (df['logFC'] > LOGFC_THRESH)  & (df['cell_type'].isin(celltype))].sort_values(['cell_type','logFC'],ascending=False)
    dn = df[(df['FDR'] < FDR_THRESH) & (df['logFC'] < -LOGFC_THRESH) & (df['cell_type'].isin(celltype))].sort_values(['cell_type','logFC'],ascending=True)
    if sortby!='LOGFC':
        up = up.sort_values(['cell_type','FDR'],ascending=True)
        dn = dn.sort_values(['cell_type','FDR'],ascending=True)
    print(f'filtering at FDR_THRESH={FDR_THRESH}, LOGFC_THRESH={LOGFC_THRESH}, length={len(up)}/{len(dn)}/{len(up)+len(dn)}')
    if print_bins:
        print(dict(zip(bins,labels)))
    return up, dn

def filterMultiLogFC(df_in,logFC_THRESH = 1, FDR_THRESH = 0.01, print_bins=None):
    import numpy as np
    import pandas as pd
    
    fixcol = ['logCPM','F','PValue','FDR']
    cols = set(df_in.columns) - set(fixcol)
    output = df_in[(np.abs(df_in[cols]).max(1) > logFC_THRESH) & (df_in['FDR'] < FDR_THRESH)]
    print(f'filtering using: logFC_THRESH={logFC_THRESH}, FDR_THRESH={FDR_THRESH}, output_size={len(output)}')

    bins=[0,1.00e-04,1.00e-03,1.00e-02,5.00e-02,1]
    labels=['****','***','**','*','ns']
    df_in['FDR_sym'] = pd.cut(df_in['FDR'], bins=bins, right=True, labels=labels)
    if print_bins:
        print(dict(zip(bins,labels)))
    return output

# plotgene 
def plotGene(gene, lum_norm_age, pro_norm_age, coltest='age_group', ax=None, lum_color='#42BCC9', pro_color='#1B758F', **kwargs):
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_style({'axes.grid' : False})
    
    if not ax: ax = plt.gca()
    sns.lineplot(data=lum_norm_age[[gene]+[coltest]], x=coltest, y=gene, errorbar='se', ax=ax, color=lum_color, marker='o', **kwargs)
    sns.lineplot(data=pro_norm_age[[gene]+[coltest]], x=coltest, y=gene, errorbar='se', ax=ax, color=pro_color, marker='o', **kwargs)
    ax.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))
    ax.tick_params(axis='both', labelsize=6)
    t = ax.yaxis.get_offset_text()
    t.set_size(6)
    ax.set_ylabel(gene, fontsize=8)
    ax.set_xlabel('', fontsize=8)
    plt.setp(ax.get_xticklabels(), rotation=45, rotation_mode="anchor", ha="right")
    plt.tight_layout()

def plotPath(pathway, mat_norm, genes, terms, ax=None, color='#1B758F'):
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns
    sns.set_style({'axes.grid' : False})

    if not ax: ax=plt.gca()
    PATH = terms[pathway]
    ss = mat_norm[set(PATH) & set(genes)].reset_index()
    ss_melt = pd.melt(ss, id_vars='age_group')
    sns.lineplot(data=ss_melt, x='age_group', y='value', errorbar='se', color=color, marker='o', ax=ax)
    # print(f'pro genes in path: {ss.shape[1]-1}')
    if not ax.get_title():
        ax.set_title(pathway);
    else:
        ax.set_title('');
    return ax

# read KEGG pathway
def read_tsvKEGG(path):
    '''read KEGG pahway from https://www.gsea-msigdb.org/gsea/msigdb/human/geneset'''
    import pandas as pd
    gl = pd.read_table(path, index_col=0)
    tmp = gl.loc['GENE_SYMBOLS'].str.extractall(r"(\w+)").values
    gl_genes = [n for sublist in tmp for n in sublist]
    return gl_genes

#  fitler genes expressed in X fraction of a cluster
def getGenesAboveRatio(adata, colname, GENE_RATIO):
    import pandas as pd
    ll = pd.DataFrame(columns=['ratio','cluster'])
    for cluster in adata.obs[colname].unique():
        idx = adata.obs[adata.obs[colname].isin([cluster])].index
        tmp = adata[idx].to_df()
        ratio = tmp.astype('bool').sum(axis=0) / tmp.shape[0]
        ratio.name = 'ratio'
        dat = pd.DataFrame(ratio)
        dat['cluster'] = cluster
        fdat = dat[dat['ratio'] > GENE_RATIO]
        print(f'{cluster} filtered genes {fdat.shape[0]}')
        ll = pd.concat([ll, fdat], axis=0)
        del tmp
    return ll.sort_values(['ratio','cluster'])

# convert edgeR output to pandas df
def edgeR2df(strvec):
    return dict(zip(strvec.names, list(strvec)))['table']

def fisher_representation(sample_size, class_in_sample, population_size, class_in_population):
    '''
    this code taken from Erick gethub
    Performs an analysis of enrichment/depletion based on observation
    in a sample. It computes a p-value given a fisher exact test.

    Parameters
    ----------
    sample_size : int
        Size of the sample obtained or number of elements
        obtained from the analysis.

    class_in_sample : int
        Number of elements of a given class that are
        contained in the sample. This is the class to be tested.

    population_size : int
        Size of the sampling space. That is, the total number
        of possible elements to be chosen when sampling.

    class_in_population : int
        Number of elements of a given class that are contained
        in the population. This is the class to be tested.

    Returns
    -------
    results : dict
        A dictionary containing the odd ratios and p-values for
        depletion and enrichment analysis.
    '''
    from scipy import stats as st
    # Computing the number of elements that are not in the same class
    nonclass_in_sample = sample_size - class_in_sample
    nonclass_in_population = population_size - class_in_population

    # Remaining elements in population after sampling
    rem_class = class_in_population - class_in_sample
    rem_nonclass = nonclass_in_population - nonclass_in_sample

    # Depletion Analysis
    depletion_odds, depletion_fisher_p_val = st.fisher_exact([[class_in_sample, rem_class],
                                                              [nonclass_in_sample, rem_nonclass]],
                                                             alternative='less')

    # Enrichment Analysis
    enrichment_odds, enrichment_fisher_p_val = st.fisher_exact([[class_in_sample, rem_class],
                                                                [nonclass_in_sample, rem_nonclass]],
                                                               alternative='greater')

    p_vals = (depletion_fisher_p_val, enrichment_fisher_p_val)
    odds = (depletion_odds, enrichment_odds)
    results = {'pval' : p_vals,
               'odds' : odds,
              }
    return results

############################
##### VISIUM functions #####
############################
def get_ax_size(ax): #ax size in px
    import matplotlib.pyplot as plt
    fig = plt.gcf()
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= fig.dpi
    height *= fig.dpi
    return width, height

# # not working 
# def addScaleBar(ax, spot_dim, scalbar_len, scalbar_wid):
#     spot_num = 65
# #     px_size = spot_num/spot_dim #um/px
#     # ss_len = scalbar_len/spot_dim #px
#     # ss_wid = scalbar_wid/spot_dim #px
#     # sizeInPx = get_ax_size(ax)
#     # relative_len = ss_len / sizeInPx[0]
#     # relative_wid = ss_wid / sizeInPx[1]
#     scalefactor = spot_dim
#     sizeInPx = get_ax_size(ax) 
#     relative_len = scalbar_len / (sizeInPx[0]*scalefactor)
#     relative_wid = scalbar_wid / (sizeInPx[1]*scalefactor)


#     x0 = 0.95 - relative_len
#     y0 = 0.03
#     x1 = relative_len
#     y1 = relative_wid
# #     print(f'{x0} {y0} {x1} {y1}')
#     axin1 = ax.inset_axes([x0, y0, x1, y1]) #can change to patch
#     axin1.set_xticks([])
#     axin1.set_yticks([])
#     for axis in ['top','bottom','left','right']:
#         axin1.spines[axis].set_linewidth(0.75)
    
#     axin1.annotate(f'{scalbar_len}', (0.5, 0.5),
#                                transform=ax.transAxes,
#                                ha='center', va='center', fontsize=6,
#                                color='black')
def cmap_map(function, cmap):
    """ Applies function (which should operate on vectors of shape 3: [r, g, b]), on colormap cmap.
    This routine will break any discontinuous points in a colormap.
    light_jet = cmap_map(lambda x: x/2 + 0.5, matplotlib.cm.jet)
    """
    import matplotlib
    import numpy as np
    import matplotlib.pyplot as plt
    cdict = cmap._segmentdata
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ('red', 'green', 'blue'):
        step_dict[key] = list(map(lambda x: x[0], cdict[key]))
    step_list = sum(step_dict.values(), [])
    step_list = np.array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step : np.array(cmap(step)[0:3])
    old_LUT = np.array(list(map(reduced_cmap, step_list)))
    new_LUT = np.array(list(map(function, old_LUT)))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i, key in enumerate(['red','green','blue']):
        this_cdict = {}
        for j, step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j, i]
            elif new_LUT[j,i] != old_LUT[j, i]:
                this_cdict[step] = new_LUT[j, i]
        colorvector = list(map(lambda x: x + (x[1], ), this_cdict.items()))
        colorvector.sort()
        cdict[key] = colorvector

    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)

def splotVis(adata, selected_clusters, figsize=(15,5), wspace=0.4,hspace=0.3, ncol=5, alpha_img=1, vmin=0, vmax='p99.5',crop=True, size=1.2, **kwargs):
    import matplotlib.pyplot as plt
    import scanpy as sc
    crop_coord = None
    
    clus2plot    = ['slide'] + selected_clusters
    spotSize     = [0] + len(selected_clusters)*[None]
    showScaleBar = [True] + len(selected_clusters)*[None]

    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(hspace=hspace, wspace=wspace)
    fig.suptitle(f"{adata.uns['metain']['donor']}_{adata.uns['metain']['section']}")

    n_cols = ncol if len(clus2plot) > ncol else len(clus2plot)
    n_rows = (len(clus2plot) // n_cols)+1 if (len(clus2plot) // n_cols)>0 else (len(clus2plot) // n_cols)
    
    if crop:
        vertices = adata.uns['metain']['capture_area_mask_vertices']
        x1 = vertices[1,1]
        x2 = vertices[3,1]
        y1 = vertices[1,0]
        y2 = y1+3850#vertices[3,0]
        crop_coord=(x1, x2, y1, y2)

    for i, (clus,spot_size, scalebar) in enumerate(zip(clus2plot, spotSize, showScaleBar)):
            ax = fig.add_subplot(n_rows, n_cols, i+1)
            f = sc.pl.spatial(adata, cmap='magma',
                              color=clus,
                              spot_size=spot_size, size=size, alpha_img=alpha_img, 
                              img_key='hires',vmin=vmin, vmax=vmax,
                              legend_loc=None, legend_fontsize=3,
                              crop_coord=crop_coord,
                            #   title=f'{clus} | {cellmaxAbundance[clus]:.1f}',
                              ax=ax, show=False, return_fig=True
                              )
            # # addScaleBar
            # if scalebar:
            #     scalbar_len = 100 #um
            #     scalbar_wid = 20
            #     spot_dim = adata.uns['spatial'][list(adata.uns['spatial'].keys())[0]]['scalefactors']['spot_diameter_fullres']
            #     addScaleBar(ax, spot_dim, scalbar_len, scalbar_wid)
    return fig

def calcRec(vertices):
    '''return: [x, y, width, height]'''
    import numpy as np
    x = vertices[0,1]
    y = vertices[0,0]
    width = vertices[1,1] - x
    height = vertices[2,0] - y
    return np.array([x, y, width, height])

# def getGenesAboveRatio(adata, colname, GENE_RATIO):
#     ll = []
#     for cluster in adata.obs[colname].unique():
#         idx = adata.obs[adata.obs[colname].isin([cluster])].index
#         tmp = adata[idx].to_df()
#         ratio = tmp.astype('bool').sum(axis=0) / tmp.shape[0]
#         ll.append(ratio[ratio > GENE_RATIO].index.tolist())
#         del tmp
#     flatlist = [l for sublist in ll for l in sublist]
#     flatlist = set(flatlist)
#     return list(flatlist)

    # # get genes expressed at least x% per cluster
    # GENE_RATIO = 0.2
    # colname = 'lineage'
    # ll = getGenesAboveRatio(adataDown, colname, GENE_RATIO)

## find low number of cells per cluster
#     # round down function
# def round_to_multiple(number, multiple):
#     return int(multiple * np.floor(number / multiple))

# target_cells={}
# MIN_CELLS = 20
# for celltype in adata.obs['celltype_asTHY']:
#     c_T21 = adata.obs[(adata.obs['celltype_asTHY'] == celltype) & (adata.obs['karyotype'] == 'T21')].index.shape[0]
#     c_2N = adata.obs[(adata.obs['celltype_asTHY'] == celltype) & (adata.obs['karyotype'] == '2n')].index.shape[0]
#     if (c_T21 > MIN_CELLS) & (c_2N > MIN_CELLS):
#         target_cells[celltype] = round_to_multiple(np.minimum(c_T21, c_2N), 10)
# print(target_cells)

# # # downsample clusters 
# # Downsample to the smallest pop size per group of T21 and 2N
# adataDown = adataDown[adataDown.obs['celltype_asTHY'].isin(target_cells.keys())]
# adataDown.var = adataDown.var[['n_counts','n_cells']]
# adataDown.obs['cell'] = adataDown.obs.index
# adatas = [adataDown[adataDown.obs['celltype_asTHY_kar'].isin([cl])] for cl in adataDown.obs['celltype_asTHY_kar'].unique()]

# for dat in adatas:
#     target = target_cells[dat.obs['celltype_asTHY'].unique()[0]]
#     if dat.n_obs > target:
#         sc.pp.subsample(dat, n_obs=target, random_state=0)
# adata_downsampled = adatas[0].concatenate(*adatas[1:])

# adataDown = adataDown[[ i in adata_downsampled.obs.cell.tolist() for i in adataDown.obs['cell'] ]]
# sc.pp.filter_genes(adataDown, min_cells=20)
# adataDown.obs['celltype_asTHY_kar'].value_counts()