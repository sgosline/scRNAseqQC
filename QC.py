import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pylab as plt

def find_inflection(ann_data, inflection_percentiles = [0,15,30,100]):
    ann_data_cumsum = np.cumsum(ann_data.obs['n_counts'])
    x_vals=np.arange(0,ann_data.obs.shape[0])
    secant_coef=ann_data_cumsum[ann_data.obs.shape[0]-1]/ann_data.obs.shape[0]
    secant_line=secant_coef*x_vals
    secant_dist=abs(ann_data_cumsum-secant_line)
    inflection_percentiles_inds = np.percentile(x_vals,inflection_percentiles).astype(int)
    inflection_points = secant_dist.argsort()[::-1]
    percentile_points = inflection_points[inflection_percentiles_inds]
    color=plt.cm.tab10(np.linspace(0,1,ann_data.obs.shape[0]))
    plt.figure(figsize=(20,10))
    plt.plot(np.array(ann_data_cumsum), label="Cumulative Sum")
    #plt.plot(np.array(secant_line), label="Secant Line")
    plt.plot(np.array(secant_dist), label="Secant Distance")
    for percentile in percentile_points:
        plt.axvline(x=percentile,ymin=0,c=color[percentile],linestyle='--',linewidth=2,label="Inflection point {}".format(percentile))
    plt.legend()
    print("Inflection point at {} for {} percentiles of greatest secant distances".format(percentile_points,inflection_percentiles))
    
def reorder_AnnData(AnnData, descending = True):
    AnnData.obs['n_counts'] = AnnData.X.sum(axis=1)
    if(descending==True):
        new_order = np.argsort(AnnData.obs['n_counts'])[::-1]
    elif(descending==False):
        new_order = np.argsort(AnnData.obs['n_counts'])[:]
    AnnData.X = AnnData.X[new_order,:].copy()
    AnnData.obs = AnnData.obs.iloc[new_order].copy()

def arcsinh_transform(AnnData, cofactor = 1000):
    AnnData.X = np.arcsinh(AnnData.X*cofactor,dtype='float')
