'''
run-qc.py
runs the scRNAseqwc tools from arguments
'''

import pandas as pd
from optparse import OptionParser
from QC import *




'''
Reads in file and returns `test_lib` object
'''
def read_file(fname):
    test_in=pd.read_csv(fname,header=0,index_col=0) #
    csr_counts=(test_in.values)
    cellIDs=test_in.index.values.astype(str)
    geneIDs=test_in.columns.values.astype(str)
    test_lib=library_data(csr_counts,cellIDs,geneIDs,sort=True) #It is preferable to sort the library sizes as downstream processes depend on this ordering.
    return test_lib



def dimReduce(test_lib,quantile=30):
    #first find inflection point
    inflections=test_lib.find_inflection()

    ##select for quantile-- is this a dictionary?
    test_dr=dimension_reduction(test_lib,inflections[quantile])

    test_dr.lib_size_normalize()
    test_dr.arcsinh_transform() #log1p is also supported by .log1p_transform
    test_dr.runPCA()
    test_dr.runUMAP() #tSNE is also supported by .runTSNE
    return test_dr




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

    lib=read_file(options['filename'])
    dimred=lib(lib,options['quantile'])



if __name__=='__main__':
    main()
