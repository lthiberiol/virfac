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
import ete2
from ete2.treeview.faces import add_face_to_node
import matplotlib.colors as colors
import matplotlib.cm as cmx
from PyQt4 import QtGui


#
# load previously processed variables
#
# dictionary of homologous groups containing which genomes are present and their respective proteins
#   groups = { 'group1': { 'genome':[ proteins ] } }
groups = load(open('homologous_groups-merged.pkl'))
# list of groups identified as part of T3SS
#   t3ss_groups = [ 'group1', 'group2', ..., 'groupN' ]
t3ss_groups = map(  lambda pos: pos.replace('&', '-'),
                    load( open( 't3ss_groups-merged.pkl' ) ) 
                 )
# list of groups described as t3ss effectors by VFDB
#   effectors = [ 'group1', 'group2', ..., 'groupN' ]
effectors = load( open( 't3ss_effector_groups.pkl' ) )
apparatus = list( set(t3ss_groups).difference( effectors ) )
# presence absence of groups in genomes DataFrame
#   pd.DataFrame(index=<genome_list>, columns=<group_list>, data=[ [0,1,0,...],...,[...] ])
presence_absence = pd.read_table('presence_absence-merged.tab', index_col=0, nrows=105)
presence_absence = presence_absence.reindex(index=tree.get_leaf_names())

sources       = pd.read_table('isolation_sources', index_col=0)
source_groups = sources.groupby('Isolation_sources')
sources       = sources['Isolation_sources']

#
# ...previously done
#
#true_positives = load( open( 'true_positives.pkl' ) )
#hmm_positives = load( open( 'hmm_positives.pkl' ) )
#merged_groups_desc = {}
#for group in t6ss_groups:
#    if group in true_positives:
#        merged_groups_desc[group] = hmm_positives[group]['desc']
#    else:
#        tmp_desc = set()
#        for tmp_group in group.split('-'):
#            if tmp_group == '14650_42938.2104' or tmp_group == '323223_93603.1872':
#                tmp_group = '14650_42938.2104-323223_93603.1872'
#            tmp_desc.add(hmm_positives[tmp_group]['desc'])
#        if len(tmp_desc) != 1:
#            print 'whaaaaat?'
#            break
#        else:
#            merged_groups_desc[group] = hmm_positives[tmp_group]['desc']
#out = open('homologous_groups_descriptions-merged.pkl', 'wb')
#dump(merged_groups_desc, out)
#out.close()

merged_groups_desc = load( open( 'homologous_groups_descriptions-merged.pkl' ) )

################################################################################
############################################################### binary heatmap #
################################################################################

apparatus_pa = presence_absence[apparatus]
effectors_pa = presence_absence[effectors]

apparatus_dist = distance.squareform( distance.pdist( apparatus_pa.T, 'jaccard' ) )
apparatus_clst = sch.linkage( apparatus_dist, method='average', metric='jaccard' )
apparatus_dend = sch.dendrogram( apparatus_clst, labels=apparatus_pa.columns )

effectors_dist = distance.squareform( distance.pdist( effectors_pa.T, 'jaccard' ) )
effectors_clst = sch.linkage( effectors_dist, method='average', metric='jaccard' )
effectors_dend = sch.dendrogram( effectors_clst, labels=effectors_pa.columns )

sorted_apparatus = []
for position in apparatus_dend[ 'leaves' ]:
    sorted_apparatus.append( apparatus_pa.columns[ position ] )
sorted_effectors = []
for position in effectors_dend[ 'leaves' ]:
    sorted_effectors.append( effectors_pa.columns[ position ] )

apparatus_pa = apparatus_pa.reindex( columns=sorted_apparatus )
effectors_pa = effectors_pa.reindex( columns=sorted_effectors )

out = open('t3ss_apparatus-presence_absence.tab', 'wb')
out.write('#Names')
out.close()
apparatus_pa.to_csv('t3ss_apparatus-presence_absence.tab', sep='\t', mode='a')
out = open('t3ss_effectors-presence_absence.tab', 'wb')
out.write('#Names')
out.close()
effectors_pa.to_csv('t3ss_effectors-presence_absence.tab', sep='\t', mode='a')

def get_color_gradient(self):
    cNorm  = colors.Normalize(vmin=0, vmax=1)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plt.get_cmap(self.colorscheme))
    color_scale = []
    for scale in np.linspace(0, 1, 201):
        hex_color = '#%02x%02x%02x' %scalarMap.to_rgba(scale)[:3]
        [r,g,b,a] = scalarMap.to_rgba(scale, bytes=True)
        color_scale.append( QtGui.QColor( r, g, b, a ) )

    return color_scale

ete2.ProfileFace.get_color_gradient = get_color_gradient

cNorm  = colors.Normalize(vmin=0, vmax=source_groups.ngroups)
scalarMap = cmx.ScalarMappable( norm=cNorm, cmap=plt.get_cmap('Paired') )
nameFaces = {}
count = 0
for isolation_source in source_groups.groups:
    nameFaces[isolation_source] = ete2.AttrFace("name", fsize=10, fgcolor='#%02x%02x%02x' %scalarMap.to_rgba(count, bytes=True)[:3])
    count +=1

profileFace = ete2.ProfileFace(1, 0.01, 0.5, len(effectors)*7.5, 14, "heatmap", colorscheme='Reds')
nameFace = ete2.AttrFace("name", fsize=10)
supportFace = ete2.AttrFace("support", fsize=10)
def myLayout(node):
    if node.is_leaf():
        add_face_to_node(nameFaces[sources[node.name]], node, 1, aligned=True)
        add_face_to_node(profileFace, node, 0, aligned=True)
        node.img_style['fgcolor'] = '#000000'
    else:
        add_face_to_node(supportFace, node, 0)
        node.img_style['size'] = 0

tree = ete2.ClusterTree('tree_puzzle/nucleotide/reconciliation/RAxML_bipartitions.reference_1000', 't3ss_effectors-presence_absence.tab', fdist=jaccard)
treeStyle                   = ete2.TreeStyle()
treeStyle.layout_fn         = myLayout
treeStyle.show_leaf_name    = False
treeStyle.draw_guiding_lines= True
treeStyle.legend.add_face(ete2.TextFace('  Isolation sources:  ', fsize=20, fgcolor='black'), column=0)
count = 1
for isolation_source in source_groups.groups:
    if count > 7:
        count = 1
    treeStyle.legend.add_face(ete2.TextFace('%s\t' %isolation_source, fsize=20, fgcolor=nameFaces[isolation_source].fgcolor), column=count)
    count += 1
treeStyle.legend_position = 3
tree_plot=tree.render('t3ss_effectors_heatmap.pdf', dpi=300, tree_style=treeStyle, units='mm')
print 'yeah'

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
    group = group.replace('&', '-')
    out.write(' <th><a href="fastas/%s.faa" target="_blank" style="color:#fff">%s<br>(%s)</a></th>' %(group[:200], merged_groups_desc[group], group))
for index in reversed(sorted_species):
    out.write('\t<tr><th style="background:#000; color:#fff; text-align:left;">%s</th>\n' %index)
    for column in sorted_genes:
        if index in groups[column]:
            tmp = []
            for protein in groups[column][index]:
                tmp.append('<a href="interproscan/%s.html">%s</a>' %(protein, protein))
            out.write('\t\t<td style="text-align:center;">%s</td>\n' %'<br>'.join(tmp))
        else:
            out.write('\t\t<td bgcolor="lightgrey"></td>\n')
    out.write('</tr>\n')
out.write('''
</table>
</body>
</html>''')
out.close()
