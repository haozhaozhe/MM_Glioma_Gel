import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import os
from glob import glob
import numpy as np
import math
from scipy import stats
from scipy.spatial import distance
import seaborn as sns

#load data
adp = sc.read('./Patch_seq_37_human_gene.h5ad')

gbm = sc.read('/gbm_4_34023.h5ad')

#plot all cells#
plot_genes = ['NTRK2', 'CACNA1E',  'KCNAB2', 
              'SHANK1', 'NEFL','GABRB1','GRIN2B',
             'ITPR1','CHRNA3','RAB11A','NRXN3','PCLO']
sc.pl.dotplot(gbm, var_names = plot_genes, groupby = 'type', swap_axes = True, cmap = 'Blues',
              standard_scale = 'var', save = '/synapes_gene_Wang_14gene.pdf')
sc.pl.dotplot(adp, var_names = plot_genes, groupby = 'Treatment', swap_axes = True,standard_scale = 'var', save = '/synapes_gene_Patch_14gene.pdf')

#plot gene score
cate_gene = ['USP8','FBXL20','RNF19A','NEDD4','FBXL20']
sc.tl.score_genes(gbm, cate_gene,  score_name='score', random_state=0, copy=False, use_raw=None)
sc.pl.dotplot(gbm, var_names = 'score',groupby = 'type', cmap = 'Blues')

sc.tl.score_genes(adp, cate_gene,  score_name='score', random_state=0, copy=False, use_raw=None)
sc.pl.dotplot(adp, var_names = 'score',groupby = 'Treatment')

cate_gene = ['ERC1','PCLO']
sc.tl.score_genes(gbm, cate_gene,  score_name='score', random_state=0, copy=False, use_raw=None)
sc.pl.dotplot(gbm, var_names = 'score',groupby = 'type', cmap = 'Blues')

sc.tl.score_genes(adp, cate_gene,  score_name='score', random_state=0, copy=False, use_raw=None)
sc.pl.dotplot(adp, var_names = 'score',groupby = 'Treatment')


#plot for each cell#
gene_list = ['MKI67','ATRX', 'EGFR', 'BRAF','NTRK2','TTYH2','CACNA1E',
             'CNTN1','HOMER1','CACNA1D','SHANK1','GRIN2B',
             ]

figsize(3,0.5)
cellID = 'A31'
adp.X = adp.layers['raw']

def plot_marker_dot(ad, cellID):
    plt.scatter(x = np.arange(0,len(gene_list)),y =[1]*len(gene_list),
                c = ad[cellID,gene_list].X, cmap = 'Reds', 
                linewidth=0.5, edgecolors = 'k', s = 80, vmax =0.005)
    plt.xticks(range(0,len(gene_list)), gene_list, rotation = 90);
    plt.ylim(0.9, 1.1)
    plt.title(cellID)
    plt.tight_layout()
    plt.savefig(f'figures/dotplot_patch_each/{cellID}_marker.png', dpi = 600,bbox_inches='tight')
    plt.savefig(f'figures/dotplot_patch_each/{cellID}_marker.pdf', dpi = 600,bbox_inches='tight')
    #ax3.set_xticklabels(gene_list, rotation = 90);
    plt.show()

for cellID in adp.obs_names:
    plot_marker_dot(adp, cellID)


#plot for tumor score
tumor_gene = ['MKI67','ATRX', 'EGFR', 'BRAF','NTRK2','TTYH2','CACNA1E',
             'CNTN1','HOMER1','CACNA1D','SHANK1','GRIN2B',
             ]
sc.tl.score_genes(adp, tumor_gene,  score_name='tumor_score', random_state=0, copy=False, use_raw=None)
sc.pl.dotplot(adp, var_names = 'tumor_score',groupby = 'Treatment', cmap = 'Blues')

#plot cross-correlation
figsize(4,3)
ephys_measure = 'EPSC_freq_Hz'
x = adp.obs.loc[:,'tumor_score']
y = adp.obs.loc[:,ephys_measure]
labels = adp.obs.loc[:,'Treatment']
df = pd.DataFrame(dict(x=x, y=y, label=labels))

x1 = adp.obs.loc[adp.obs['Treatment']=='Tumor','tumor_score']
y1 = adp.obs.loc[adp.obs['Treatment']=='Tumor',ephys_measure]
z1 = numpy.polyfit(x1, y1, 1,)
p1 = numpy.poly1d(z1)


rho_humk, pval_humk = stats.spearmanr(x1, y1)

groups = df.groupby('label')

# Plot
fig, ax = plt.subplots()
ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
for name, group in groups:
    ax.plot(group.x, group.y, marker='o', linestyle='', ms=8, label=name)

ax.text (0.15,1500,'r = '+str(format(rho_humk,'.3f')), color = 'black')
ax.text (0.15,1300,'p = '+str(format(pval_humk,'.1e')), color = 'black')

ax.legend()
ax.plot(x1,p1(x1),"-",color = "orange")
plt.title(ephys_measure)
plt.savefig(f'./{ephys_measure}.pdf')


# In[77]:


figsize(4,3)
ephys_measure = 'baseline_pA'

x = adp.obs.loc[:,'tumor_score']
y = adp.obs.loc[:,ephys_measure]
labels = adp.obs.loc[:,'Treatment']
df = pd.DataFrame(dict(x=x, y=y, label=labels))

x1 = adp.obs.loc[adp.obs['Treatment']=='Tumor','tumor_score']
y1 = adp.obs.loc[adp.obs['Treatment']=='Tumor',ephys_measure]
z1 = numpy.polyfit(x1, y1, 1,)
p1 = numpy.poly1d(z1)


rho_humk, pval_humk = stats.spearmanr(x1, y1)

groups = df.groupby('label')

# Plot
fig, ax = plt.subplots()
ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
for name, group in groups:
    ax.plot(group.x, group.y, marker='o', linestyle='', ms=8, label=name)

ax.text (0.15,-200,'r = '+str(format(rho_humk,'.3f')), color = 'black')
ax.text (0.15,-220,'p = '+str(format(pval_humk,'.1e')), color = 'black')

ax.legend()
ax.plot(x1,p1(x1),"-",color = "orange")
plt.title(ephys_measure)
plt.savefig(f'./{ephys_measure}.pdf')


figsize(4,3)
ephys_measure = 'EPSC_Amp_pA'

x = adp.obs.loc[:,'tumor_score']
y = adp.obs.loc[:,ephys_measure]
labels = adp.obs.loc[:,'Treatment']
df = pd.DataFrame(dict(x=x, y=y, label=labels))

x1 = adp.obs.loc[adp.obs['Treatment']=='Tumor','tumor_score']
y1 = adp.obs.loc[adp.obs['Treatment']=='Tumor',ephys_measure]
z1 = numpy.polyfit(x1, y1, 1,)
p1 = numpy.poly1d(z1)


rho_humk, pval_humk = stats.spearmanr(x1, y1)

groups = df.groupby('label')


fig, ax = plt.subplots()
ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
for name, group in groups:
    ax.plot(group.x, group.y, marker='o', linestyle='', ms=8, label=name)

ax.text (0.15,-50,'r = '+str(format(rho_humk,'.3f')), color = 'black')
ax.text (0.15,-55,'p = '+str(format(pval_humk,'.1e')), color = 'black')

ax.legend()
ax.plot(x1,p1(x1),"-",color = "orange")
plt.title(ephys_measure)
plt.savefig(f'./{ephys_measure}.pdf')



