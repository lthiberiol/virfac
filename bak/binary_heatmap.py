#!/usr/bin/python
#coding: utf-8

from matplotlib import use
use('Agg')
import pandas as pd
import numpy as np
import pylab as pl
import matplotlib.gridspec as gridspec
import scipy.spatial.distance as distance
import scipy.cluster.hierarchy as sch

binary_df = pd.read_table('presence_absence.tab', index_col=0)
to_remove = []
for group in binary_df.columns:
    if binary_df[group].sum() > len(binary_df.index)*0.75:
        to_remove.append(group)
binary_df = binary_df.drop(to_remove, axis=1)

sch.set_link_color_palette(['black'])

species_dist = distance.squareform(distance.pdist(binary_df, 'jaccard'))
species_clst = sch.linkage(species_dist, method='average')
species_dend = sch.dendrogram(species_clst, labels=binary_df.index)

genes_dist = distance.squareform(distance.pdist(binary_df.T, 'jaccard'))
genes_clst = sch.linkage(genes_dist, method='average')
genes_dend = sch.dendrogram(genes_clst, labels=binary_df.columns)

sorted_species = []
for position in species_dend['leaves']:
    sorted_species.append(binary_df.index[position])
sorted_genes = []
for position in genes_dend['leaves']:
    sorted_genes.append(binary_df.columns[position])

matrix = []
for index in sorted_species:
    matrix.append([])
    for column in sorted_genes:
        matrix[-1].append(binary_df.loc[index][column])
half_matrix = []
for line in matrix:
    half_matrix.append(map(lambda x: x/2.0, line))

fig = pl.figure(figsize=(40,20))
gs = gridspec.GridSpec(1,2,wspace=0.0,width_ratios=[0.25, 1])

dend_ax = fig.add_subplot(gs[0,0])
dend    = sch.dendrogram(species_clst, color_threshold=np.inf, orientation='right')
dend_ax.set_xticks([])
for sp in dend_ax.spines.values():
    sp.set_visible(False)

heatmap_ax = fig.add_subplot(gs[0,1])
heatmap    = heatmap_ax.imshow(matrix, interpolation='nearest', aspect='auto', origin='lower', cmap='Blues')

heatmap_ax.set_yticks(np.arange(binary_df.shape[0]))
heatmap_ax.yaxis.set_ticks_position('right')
heatmap_ax.set_yticklabels(sorted_species)
heatmap_ax.set_xticks([])
heatmap_ax.grid(axis='y')

fig.tight_layout()
fig.savefig('ipad_sync/genomes_hg-heatmap-complete.pdf')
