__author__ = 'Thiberio'
from matplotlib import use
use('Agg')
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.spatial.distance as distance
import scipy.cluster.hierarchy as sch
import pandas as pd
import numpy as np
from os import chdir, listdir
from pickle import load, dump

chdir('/Volumes/Macintosh HD 2/thiberio/virulence_factors')

groups = load(open('homologous_groups-merged.pkl'))
t3ss_groups = load(open('t3ss_groups-merged.pkl'))
df = pd.read_table('presence_absence-merged.tab', index_col=0, nrows=105)
df = df[t3ss_groups]

merged_groups_desc = {}
for group in t3ss_groups:
    if '&' not in group:
        merged_groups_desc[group] = groups_desc[group]
    else:
        tmp_desc = set()
        for tmp_group in group.split('&'):
            tmp_desc.add(groups_desc[tmp_group])
        if len(tmp_desc) != 1:
            print 'whaaaaat?'
            break
        else:
            merged_groups_desc[group] = groups_desc[tmp_group]
out = open('homologous_groups_descriptions-merged.pkl', 'wb')
dump(merged_groups_desc, out)
out.close()

################################################################################
############################################################### binary heatmap #
################################################################################
species_dist = distance.squareform(distance.pdist(df, 'jaccard'))
species_clst = sch.linkage(species_dist, method='average')
species_dend = sch.dendrogram(species_clst, labels=df.index)

genes_dist = distance.squareform(distance.pdist(df.T, 'jaccard'))
genes_clst = sch.linkage(genes_dist, method='average')
genes_dend = sch.dendrogram(genes_clst, labels=df.columns)

sorted_species = []
for position in species_dend['leaves']:
    sorted_species.append(df.index[position])
sorted_genes = []
for position in genes_dend['leaves']:
    sorted_genes.append(df.columns[position])

matrix = []
for index in sorted_species:
    matrix.append([])
    for column in sorted_genes:
        matrix[-1].append(df.loc[index][column])
half_matrix = []
for line in matrix:
    half_matrix.append(map(lambda x: x/2.0, line))

fig = plt.figure(figsize=(40,20), dpi=300)
gs = gridspec.GridSpec(1,2,wspace=0.0,width_ratios=[0.25, 1])
sch.set_link_color_palette(['black'])

dend_ax = fig.add_subplot(gs[0,0])
dend    = sch.dendrogram(species_clst, color_threshold=np.inf, orientation='right')
dend_ax.set_xticks([])
for sp in dend_ax.spines.values():
    sp.set_visible(False)

heatmap_ax = fig.add_subplot(gs[0,1])
heatmap    = heatmap_ax.imshow(matrix, interpolation='nearest', aspect='auto', origin='lower', cmap='Blues')

heatmap_ax.set_yticks(np.arange(df.shape[0]))
heatmap_ax.yaxis.set_ticks_position('right')
heatmap_ax.set_yticklabels(sorted_species)
heatmap_ax.set_xticks([])
heatmap_ax.grid(axis='y')

fig.tight_layout()
fig.savefig('ipad_sync/genomes_hg-heatmap-average-merged.pdf')


################################################################################
################################################################# HTML version #
################################################################################
out = open('presence_absence.html', 'wb')
out.write('''
<html>
<body>
<table border="1">
<tr style="background:#000; color:#fff; "> <th></th>\n''' )
for group in sorted_genes:
    out.write(' <th><a href="homologous_groups/%s.faa" target="_blank" style="color:#fff">%s</a></th>' %(group, merged_groups_desc[group]))
for index in reversed(sorted_species):
    out.write('\t<tr><th style="background:#000; color:#fff; text-align:left;">%s</th>\n' %index)
    for column in sorted_genes:
        if index in groups[column]:
            tmp = []
            for protein in groups[column][index]:
                tmp.append('<a href="interproscan/%s|%s.html">%s</a>' %(protein, index, protein))
            out.write('\t\t<td style="text-align:center;">%s</td>\n' %'<br>'.join(tmp))
        else:
            out.write('\t\t<td bgcolor="lightgrey"></td>\n')
    out.write('</tr>\n')
out.write('''
</table>
</body>
</html>''')
out.close()