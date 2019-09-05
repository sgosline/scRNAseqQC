import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pylab as plt
from scipy import stats

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
    #save to file
    if(output_prefix!=''):
        plt.savefig(output_prefix+'_inflectionCheck.png',bbox_inches='tight')
    else:
        plt.show()
    print("Inflection point at {} for {} percentiles of greatest secant distances".format(percentile_points,inflection_percentiles))
    #SJCG: added the ability to return a dictionary of points
    return(dict(zip(inflection_percentiles, percentile_points)))
    
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

def cluster_summary_stats(AnnData,raw=False):
    cluster_means = np.zeros((len(np.unique(AnnData.obs['louvain'])),AnnData.n_vars))
    cluster_medians = np.zeros((len(np.unique(AnnData.obs['louvain'])),AnnData.n_vars))
    cluster_stdev = np.zeros((len(np.unique(AnnData.obs['louvain'])),AnnData.n_vars))
    if(raw == True):
        for cluster in range(len(np.unique(AnnData.obs['louvain']))):
            cluster_means[cluster]=np.array(np.mean(AnnData[AnnData.obs['louvain'].isin([str(cluster)])].X,axis = 0))
            cluster_medians[cluster]=np.array(np.median(AnnData[AnnData.obs['louvain'].isin([str(cluster)])].X,axis = 0))
            cluster_stdev[cluster]=np.array(np.std(AnnData[AnnData.obs['louvain'].isin([str(cluster)])].X,axis = 0))
    elif(raw == False):    
        for cluster in range(len(np.unique(AnnData.obs['louvain']))):
            cluster_means[cluster]=np.array(np.mean(AnnData[AnnData.obs['louvain'].isin([str(cluster)])].raw.X,axis = 0))
            cluster_medians[cluster]=np.array(np.median(AnnData[AnnData.obs['louvain'].isin([str(cluster)])].raw.X,axis = 0))
            cluster_stdev[cluster]=np.array(np.std(AnnData[AnnData.obs['louvain'].isin([str(cluster)])].raw.X,axis = 0))
    AnnData.layers['Cluster_Medians'] = np.array(cluster_medians[AnnData.obs['louvain'].astype(int)])
    AnnData.layers['Cluster_Means'] = cluster_means[AnnData.obs['louvain'].astype(int)]
    AnnData.layers['Cluster_Stdevs'] = cluster_stdev[AnnData.obs['louvain'].astype(int)]

def cluster_wilcoxon_rank_sum(AnnData,feature_list,alternative='greater'):
    cluster_list = AnnData.obs['louvain']
    p_values = np.zeros((len(np.unique(cluster_list)),len(feature_list)))
    for cluster in range(len(np.unique(cluster_list))):
        for feature in range(len(feature_list)):
            p_values[cluster,feature]=stats.mannwhitneyu(AnnData[cluster_list.isin([str(cluster)])].obs_vector(feature_list[feature]),AnnData.obs_vector(feature_list[feature]),alternative='greater',use_continuity=True)[1]
    AnnData.uns['Cluster_p_values'] = pd.DataFrame(p_values,np.arange(len(np.unique(cluster_list))),feature_list)

def cluster_p_threshold(AnnData,threshold = 0.05):
    for columns in AnnData.uns['Cluster_p_values']:
        AnnData.obs[columns+'_threshold'] = (AnnData.uns['Cluster_p_values']<threshold)[columns][AnnData.obs['louvain'].astype(int)].astype(int).values
        AnnData.obs[columns+'_enrichment'] = (AnnData.uns['Cluster_p_values'])[columns][AnnData.obs['louvain'].astype(int)].values