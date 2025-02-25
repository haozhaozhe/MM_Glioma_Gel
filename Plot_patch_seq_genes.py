#!/usr/bin/env python
# coding: utf-8

# In[11]:


get_ipython().run_line_magic('pylab', 'inline')

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


# # read data

# In[2]:


adp = sc.read('./Patch_seq_37_human_gene.h5ad')


# In[3]:


gbm = sc.read('/gbm_4_34023.h5ad')


# # group genes

# In[13]:


plot_genes = ['NTRK2', 'CACNA1E',  'KCNAB2', 
              'SHANK1', 'NEFL','GABRB1','GRIN2B',
             'ITPR1','CHRNA3','RAB11A','NRXN3','PCLO']
sc.pl.dotplot(gbm, var_names = plot_genes, groupby = 'type', swap_axes = True, cmap = 'Blues',
              standard_scale = 'var', save = '/synapes_gene_Wang_14gene.pdf')
sc.pl.dotplot(adp, var_names = plot_genes, groupby = 'Treatment', swap_axes = True,standard_scale = 'var', save = '/synapes_gene_Patch_14gene.pdf')


# In[15]:


cate_gene = ['USP8','FBXL20','RNF19A','NEDD4','FBXL20']
sc.tl.score_genes(gbm, cate_gene,  score_name='score', random_state=0, copy=False, use_raw=None)
sc.pl.dotplot(gbm, var_names = 'score',groupby = 'type', cmap = 'Blues')

sc.tl.score_genes(adp, cate_gene,  score_name='score', random_state=0, copy=False, use_raw=None)
sc.pl.dotplot(adp, var_names = 'score',groupby = 'Treatment')


# In[16]:


cate_gene = ['ERC1','PCLO']
sc.tl.score_genes(gbm, cate_gene,  score_name='score', random_state=0, copy=False, use_raw=None)
sc.pl.dotplot(gbm, var_names = 'score',groupby = 'type', cmap = 'Blues')

sc.tl.score_genes(adp, cate_gene,  score_name='score', random_state=0, copy=False, use_raw=None)
sc.pl.dotplot(adp, var_names = 'score',groupby = 'Treatment')


# # each cell

# In[18]:


gene_list = ['NES','PSD95','CHI3L1','PROM1','A2B5','FUT4','L1CAM','ATRX', 'TERT', 'H3', 'EGFR', 'BRAF','MKI67','SOX2','HOMER2','HOMER1','SHANK1','SHANK2','SHANK3',
'TTYH2','CNTN1','CNTN2','CNTN3','NTRK2','NTRK3','GJA1','CACNA1A','CACNA1B','CACNA1C','CACNA1D',
'CACNA1E','CACNA1F','CACNA1G','CACNA1I','CACNA2D1','CACNA2D2','CACNA2D3','CACNB1','CACNB2',
'CACNB3','CACNB4','CACNG2','CACNG3','CACNG8','GRIA1','GRIA2','GRIA3','GRIA4','GRIN1','GRIK1',
'GRIN2A','GRIN2B','GRIN2C','GRM1','GRM2','GRM5','GRM7','GRIN3A',]
#gene_list =  [s.capitalize() for s in gene_list]
#checked, these genes are the same between mouse and human

gene_list = ['MKI67','ATRX', 'EGFR', 'BRAF','NTRK2','TTYH2','CACNA1E',
             'CNTN1','HOMER1','CACNA1D','SHANK1','GRIN2B',
             ]
#gene_list =  [s.capitalize() for s in gene_list]
#gene_list = list(set(gene_list)&set(adp.var_names))

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


# # score correlation

# In[65]:


tumor_gene = ['MKI67','ATRX', 'EGFR', 'BRAF','NTRK2','TTYH2','CACNA1E',
             'CNTN1','HOMER1','CACNA1D','SHANK1','GRIN2B',
             ]
sc.tl.score_genes(adp, tumor_gene,  score_name='tumor_score', random_state=0, copy=False, use_raw=None)
sc.pl.dotplot(adp, var_names = 'tumor_score',groupby = 'Treatment', cmap = 'Blues')


# In[74]:


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


# In[80]:


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

# Plot
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


# In[ ]:




