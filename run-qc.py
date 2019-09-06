'''
run-qc.py
runs the scRNAseqwc tools from arguments
'''

import pandas as pd
from optparse import OptionParser
from QC import *




'''
Reads in file and returns `test_lib` object
now in anndata form
'''
def read_file(fname):
   # test_in=pd.read_csv(fname,header=0,index_col=0) #
   # csr_counts=(test_in.values)
   # cellIDs=test_in.index.values.astype(str)
   # geneIDs=test_in.columns.values.astype(str)
   # test_lib=library_data(csr_counts,cellIDs,geneIDs,sort=True) #It is preferable to sort the library sizes as downstream processes depend on this ordering.
    test_lib=sc.read_csv(fname)
    reorder_AnnData(test_lib) #quick reordering of AnnData object, since our inflection point analysis assumes that libraries are sorted from highest to lowest quality
    adata.raw = adata #checkpoint before normalizing and scaling data

    return test_lib

def calcGeneFracs(adata):
    #create new "observation" as percent mito and mean mito
    mito_genes = adata.var_names.str.startswith('mt-')
    adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1) / adata.obs['n_counts']
    adata.obs['mean_mito'] = np.mean(adata[:, mito_genes].X, axis=1)

    #create new "observation" as percent hematopoietic and mean hematopoietic
    hematopoietic_genes = adata.var_names.str.contains('Hbb-bs$|Cd14$|Jchain$|NAGK$')
    adata.obs['percent_hematopoietic'] = np.sum(adata[:, hematopoietic_genes].X, axis=1) / adata.obs['n_counts']
    adata.obs['mean_hematopoietic'] = np.mean(adata[:, hematopoietic_genes].X, axis=1)

    #It's important to examine multiple methods of compiling genes into a single meta-feature, as percents may be biased towards the number of genes used in the summation while means may be thrown off by outliers.
    return adata


'''
identify point of cell by gene inflection and remove low-read cells
'''
def inflectionGate(test_lib,quantile=30,outputprefix='test'):
    #first find inflection point
    inflections=test_lib.find_inflection(output_prefix=outputprefix)

    ##select for quantile-- is this a dictionary?
#    test_dr=dimension_reduction(test_lib,inflections[quantile])
    #updated to use scanpy
    test_dr = sc.pp.filter_cells(test_lib,mincounts=test_lib.obs.iloc[infelctions[quantile]].ncounts)
    test_dr.obs['ranked_n_counts'] = np.argsort(adata.obs['n_counts'])
    return test_dr

'''
normalize and reduce dimensions via umap
'''
def normalizeAndDimReduce(test_dr):
    sc.pp.normalize_total(test_dr) #each gene count value is divided by the total number of counts for that respective cell
    sc.pp.log1p(test_dr) #log1p normalization
    #same as above but wrapped up in one function and utilizes the arcsinh transform instead of log1p, commented out for convenience
    #arcsinh_transform(adata)
    sc.pp.scale(test_dr) #scaling by variance and centering to zero for visualization
    sc.tl.pca(test_dr) #performing PCA
    sc.pl.pca(test_dr, color=['n_counts','percent_mito'],components=['1,2','1,3','2,3'],ncols=3)
    sc.pl.pca_variance_ratio(test_dr, log=True) #examine amount of variance captured by each principal component
    sc.pp.neighbors(test_dr,n_neighbors=30, n_pcs=20) #UMAP requires this neighborhood calculation first, will give deprecation warnings
    sc.tl.umap(test_dr) #perform UMAP
    sc.pl.umap(test_dr,color=['Krt20','Muc2','Cd47','Rpl27a','Jchain','Hbb-bs','Cd14','NAGK','percent_mito','n_counts']) #plot marker genes to determine which clusters to keep or gate out
    return test_dr

'''
cluster the data
'''
def cluster(adata,res=1.6):
    sc.tl.louvain(adata,resolution=res)
    #the resolution value here determines the granularity of clustering
    #higher resolution = smaller, refined clusters
    #lower resolution = larger, coarse grained clusters
    sc.pl.umap(adata, color=['louvain'],wspace=0.5,palette='tab20')
    return adata

def doGating(test_dr,umapOrTsne='UMAP',density=10,delta=0.5):
    test_gate = gate_visualize(test_dr)
    if(umapOrTsne=='UMAP'):
        test_gate.runDPC(test_gate.UMAP,density,delta) #first value = density cutoff, second value = delta/a.u. cutoff
    else:
        test_gate.runDPC(test_gate.TSNE,density,delta)
    return test_gate



'''
main
'''
def main():
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename", help="name of CSV to -cells are rows, genes are columns")
    parser.add_option("-q", "--quantile", dest="quant",
                  help="Quantile to select for inflection point")
    parser.add_option('-o','--output',dest='output',help='prefix of output file')


    (options, args) = parser.parse_args()

     #read in file
    lib = read_file(options['filename'])
    #calculate genes of interest
    lib = calcGeneFracs(lib)

    #do the inflection point testing and the clustering
    dimred=dimReduce(lib,options['quantile'],outputprefix=options['output'])



if __name__=='__main__':
    main()
