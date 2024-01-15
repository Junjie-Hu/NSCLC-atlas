import numpy as np
import scanpy as sc
import scvelo as scvh
import tensorflow as tf
import pandas as pd
import unitvelo as utv
import phylovelo as pv
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

print(tf.config.list_physical_devices('GPU'))

adata = sc.read("GB_CAF_velo.h5ad")

scv.pp.filter_and_normalize(adata,min_shared_counts=10)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
adata.var['highly_variable'] = True
velo_config = utv.config.Configuration()
velo_config.R2_ADJUST = True
velo_config.IROOT = None
velo_config.FIT_OPTION = '1'
velo_config.AGENES_R2 = 1
adata = utv.run_model(adata, 'Subtype', config_file=velo_config)

scv.tl.latent_time(adata,min_likelihood=None)

# plot
adata.uns['Subtype_colors'] = ['#B3DE69','#FB8072','#FCCDE5']
scv.pl.velocity_embedding_stream(adata, basis='X_umap',color='Subtype',
                                 s=255,alpha=0.7,density=1,arrow_size=1.5,smooth=1,
                                 legend_loc='right margin',legend_fontsize=5,
                                 save="CAF_pre_unitveloumap_stream_new.pdf", figsize=(6.5,5)
                                )

celltype = ['ADH1B+ CAF', 'FAP+aSMA+ CAF', 'FAP+aSMA- CAF']
cmap = ['#B3DE69','#FB8072','#FCCDE5']
color_map = dict(zip(['ADH1B+ CAF', 'FAP+aSMA+ CAF', 'FAP+aSMA- CAF'], [0, 1, 2]))
state_map = dict(zip(['ADH1B+ CAF', 'FAP+aSMA+ CAF', 'FAP+aSMA- CAF'], ['ADH1B+ CAF', 'FAP+aSMA+ CAF', 'FAP+aSMA- CAF']))

name = 'GB_CAF_unitvelo'
phytime_bar = dict()
for i in state_map:
    phytime_bar[state_map[i]] = adata.obs.latent_time.to_numpy()[np.where(np.array(adata.obs['Subtype'])==str(i))]
hist_data = []
hist_labels = []
hist_colors = []
for i in celltype:
    hist_data.append(phytime_bar[state_map[i]])
    hist_labels.append(state_map[i])
    hist_colors.append(cmap[color_map[i]])

hd = plt.hist(hist_data, label=hist_labels, color=hist_colors,alpha=0.7)
fig, ax = plt.subplots(figsize=(8,6))

pv.ana_utils.mullerplot(hd[0],hist_labels, hist_colors, absolute=0,alpha=0.6, ax=ax)
ax.set_xlim(0, 8.5)
ax.set_ylim(0, 1)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('Latent time', font='arial', fontsize=15)
ax.set_ylabel('Fraction', fontsize=15)
plt.savefig('/figures/'+name+'_muller.png', format='png')
