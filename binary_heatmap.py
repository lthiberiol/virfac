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
import ete3
from ete3.treeview.faces import add_face_to_node
import matplotlib.colors as colors
import matplotlib.cm as cmx
from PyQt4 import QtGui
import seaborn as sns
from itertools import combinations
import multiprocessing
from scipy.stats import pearsonr

#
# load previously processed variables
#
# dictionary of homologous groups containing which genomes are present and their respective proteins
#   groups = { 'group1': { 'genome':[ proteins ] } }
groups = load(open('homologous_groups-merged.pkl'))

# list of groups identified as part of T3SS
descriptions = {}
for group, desc in load( open('homologous_groups_descriptions-merged.pkl') ).items():
        descriptions[ group.replace('&', '-') ] = desc 

new_descriptions = [ ['70959_26797.2785' ,'axoU'],      ['252396_642.133.peg.1627-245405_42950.3198' ,'aopS'], 
                     ['70378_26797.2204' ,'aopT'],      ['70973_26797.2799' ,'aseG'], 
                     ['321199_93602.2926' ,'aopX'],     ['722_642.148.peg.724' ,'aopQ'], 
                     ['320898_93602.956' ,'ateA'],      ['321039_93602.2082-165230_103.4199' ,'ateB'], 
                     ['1896_642.148.peg.1898' ,'aopH2'],['242127_46596.4542','None'],
                     ['321038_93602.2081', 'atcB'],     ['274108_35496.2207', 'tccC3'],
                     ['198566_18299.3085', 'tccC3'] ]

for group, desc in new_descriptions:
    descriptions[group] = desc
                                                                                                    
revised_effector_descriptions = set( ['aexu', 'aopx', 'ateb', 'atea', 'aext', 'aoph2', 'aoph', 'axou', 'ati2', 'tccc3', 'aopo', 'aseg', 'aops', 'aopq', 'aopt', 'xopac', 'spvb', 'exoy', 'ospf', 'exou'] )
true_positive_effectors  = []
for group in descriptions.keys():
    if descriptions[group].lower() in revised_effector_descriptions:
        true_positive_effectors.append(group)

pa = pd.read_table( 'presence_absence-merged.tab', index_col=0 )
presence_absence = pd.read_table( 'presence_absence-merged.tab', index_col=0 )
presence_absence = presence_absence.reindex( columns=true_positive_effectors )

description_pa                = pd.DataFrame( index=presence_absence.index, columns=revised_effector_descriptions, data=0)
for group in true_positive_effectors:
    present_in_genomes = presence_absence[presence_absence[group] == 1].index
    description_pa.loc[present_in_genomes, descriptions[group].lower()] = [1] * len(present_in_genomes)

aopP_pa = []
for genome in description_pa.index:
    if genome in ['A_salmonicida_01B526', 'A_salmonicida_smithia_CIP104757', 'A_salmonicida_masoucida_CIP103210', 'A_salmonicida_achromogenes_CIP104001', 'A_salmonicida_CIP103209T', 'A_salmonicida_AS03']:
        aopP_pa.append(1)
    else:
        aopP_pa.append(0)
description_pa.loc[:, 'aopp'] = pd.Series( data=aopP_pa, index=description_pa.index)

sources       = pd.read_table('isolation_sources', index_col=0)
source_groups = sources.groupby('Isolation_sources')
sources       = sources['Isolation_sources']

x = []
y = []
for source in source_groups.groups.keys():
    tmp = description_pa.loc[ source_groups.get_group(source).index ].T.sum()
    x.append( source )
    y.append( tmp.values )
    print source, tmp.mean(), tmp.std()

fig, ax = plt.subplots()
sns.violinplot( vals=y, names=x, ax=ax )
fig.tight_layout()
fig.savefig('yeah.pdf')

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

def pearson_wrapper((x,y)):
    return list(pearsonr(x,y))
print 'yeah'

pearson_wrapper( description_pa.T.values[:5] )
pool = multiprocessing.Pool(processes=10)
pearson_output = pool.map(pearson_wrapper, list( combinations(description_pa.T.values, 2)) )
pool.close()
pool.join()
print 'yeah'

pearson_output = np.asarray(pearson_output)

r = distance.squareform(pearson_output[:, 0])
p = distance.squareform(pearson_output[:, 1])
mask = p > 0.05

fig, ax = plt.subplots()
sns.heatmap(r, ax=ax, cmap='RdBu', xticklabels=sorted_effectors, yticklabels=sorted_effectors)
fig.tight_layout()
fig.savefig('eff_corr_heatmap3.png', dpi=300, figsize=[10,10])


effectors_dist = distance.squareform( distance.pdist( description_pa.T, 'correlation' ) )
effectors_clst = sch.linkage( effectors_dist, method='average', metric='correlation' )
effectors_dend = sch.dendrogram( effectors_clst, labels=description_pa.columns, orientation='left' )

sorted_effectors = []
for position in effectors_dend[ 'leaves' ]:
    sorted_effectors.append( description_pa.columns[ position ] )

description_pa = description_pa.reindex( columns=sorted_effectors )

out = open('t3ss_effectors_descriptions-presence_absence.tab2', 'wb')
out.write('#Names')
out.close()
description_pa.to_csv('t3ss_effectors_descriptions-presence_absence.tab2', sep='\t', mode='a')

def get_color_gradient(self):
    from PyQt4 import QtGui
    import matplotlib.colors as colors
    import matplotlib.cm as cmx

    cNorm  = colors.Normalize(vmin=0, vmax=1)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plt.get_cmap(self.colorscheme))
    color_scale = []
    for scale in np.linspace(0, 1, 201):
        hex_color = '#%02x%02x%02x' %scalarMap.to_rgba(scale)[:3]
        [r,g,b,a] = scalarMap.to_rgba(scale, bytes=True)
        color_scale.append( QtGui.QColor( r, g, b, a ) )

    return color_scale

ete3.ProfileFace.get_color_gradient = get_color_gradient

cNorm  = colors.Normalize(vmin=0, vmax=source_groups.ngroups)
scalarMap = cmx.ScalarMappable( norm=cNorm, cmap=plt.get_cmap('Paired') )
nameFaces = {}
sourceFaces = {}
count = 0
for isolation_source in source_groups.groups:
    sourceFaces[isolation_source] = ete3.CircleFace(6, '#%02x%02x%02x' %scalarMap.to_rgba(count, bytes=True)[:3])
    nameFaces[isolation_source] = ete2.AttrFace("name", fsize=10, fgcolor='#%02x%02x%02x' %scalarMap.to_rgba(count, bytes=True)[:3])
    count +=1

#
# without heatmap
#
nameFace = ete2.AttrFace("name", fsize=10)
nameFace_no_isolation_source = ete2.AttrFace("name", fsize=10, fgcolor='#000000' )
supportFace = ete2.AttrFace("support", fsize=10)
def myLayout(node):
    if node.is_leaf():
        add_face_to_node(nameFaces[sources[node.name]], node, 1, aligned=True)
        #add_face_to_node(nameFace_no_isolation_source, node, 1, aligned=True)
        node.name = node.name.replace('_', ' ')
        node.img_style['size'] = 0
    elif node.is_root():
        node.img_style['size'] = 0
    else:
        if node.support >= 95:
            node.img_style['size'] = 7
            node.img_style['fgcolor'] = '#791205'
        else:
            node.img_style['size'] = 0

tree = ete2.Tree( 'tree_puzzle/nucleotide/concatenation/RAxML_fastTreeSH_Support.sh_like-unrooted' )
aeromonas_outgroup = tree.get_common_ancestor('A_simiae_CIP107798T', 'A_schubertii_CECT4240T', 'A_diversa_CECT4254T')
tree.set_outgroup( aeromonas_outgroup )
treeStyle                   = ete2.TreeStyle()
treeStyle.layout_fn         = myLayout
treeStyle.show_leaf_name    = False
treeStyle.draw_guiding_lines= True
treeStyle.legend.add_face(ete2.TextFace('  Isolation sources:  ', fsize=20, fgcolor='black'), column=0)
count = 1
for isolation_source in source_groups.groups:
    if count > 7:
        count = 1
    treeStyle.legend.add_face(ete2.TextFace('%s    ' %isolation_source, fsize=20, fgcolor=nameFaces[isolation_source].fgcolor), column=count)
    count += 1
treeStyle.legend_position = 3
tree_plot=tree.render('cluster0_concatenation-isolation_sources.pdf', dpi=300, tree_style=treeStyle, units='mm')

#
# with heatmap
#
heatmap = ete3.ProfileFace(1, 0.01, 0.5, description_pa.shape[1]*15, 14, "heatmap", colorscheme='Greys')
nameFace = ete3.AttrFace("name", fsize=10)
supportFace = ete3.AttrFace("support", fsize=10)
def myLayout(node):
    if node.is_leaf():
        add_face_to_node(nameFace, node, 2, aligned=True)
        add_face_to_node(heatmap, node, 0, aligned=True)
        add_face_to_node(sourceFaces[sources[node.name]], node, 1, aligned=True)
        node.name = node.name.replace('_', ' ')
        node.img_style['size'] = 0
    elif node.is_root():
        node.img_style['size'] = 0
    else:
        if node.support >= 95:
            node.img_style['size'] = 7
            node.img_style['fgcolor'] = '#791205'
        else:
            node.img_style['size'] = 0
#def myLayout(node):
#    node.img_style['size'] = 0
#    if node.is_leaf():
#        add_face_to_node(nameFaces[sources[node.name]], node, 1, aligned=True)
#        add_face_to_node(profileFace, node, 0, aligned=True)
#        node.name = node.name.replace('_', ' ')

tree = ete3.ClusterTree( 'tree_puzzle/nucleotide/concatenation/RAxML_fastTreeSH_Support.sh_like-unrooted', 't3ss_effectors_descriptions-presence_absence.tab2')
aeromonas_outgroup = tree.get_common_ancestor('A_simiae_CIP107798T', 'A_schubertii_CECT4240T', 'A_diversa_CECT4254T')
tree.set_outgroup( aeromonas_outgroup )
treeStyle                   = ete3.TreeStyle()
treeStyle.layout_fn         = myLayout
treeStyle.show_leaf_name    = False
treeStyle.draw_guiding_lines= True
#treeStyle.legend.add_face(ete3.TextFace('  Isolation sources:  ', fsize=20, fgcolor='black'), column=0)
#count = 1
#for isolation_source in source_groups.groups:
#    if count > 7:
#        count = 1
#    treeStyle.legend.add_face(ete3.TextFace('%s    ' %isolation_source, fsize=20, fgcolor=nameFaces[isolation_source].fgcolor), column=count)
#    count += 1
#treeStyle.legend_position = 3
tree_plot=tree.render('cluster0_concatenation-isolation_sources-desc_heatmap2.svg', dpi=800, tree_style=treeStyle, units='mm')
print 'yeah'

################################################################################
############################## evolutionary distances within isolation sources #
################################################################################
reference_dist = pd.read_table('tree_puzzle/nucleotide/distances/320179_93602.2343.phy.dist.tab', index_col=0)

new_index = []
for index in reference_dist.index:
    new_index.append( index.split('|')[1] )

reference_dist.index   = new_index
reference_dist.columns = new_index

sns.set(style="ticks")
fig, ax = plt.subplots( dpi=300 )
for isolation_source in source_groups.groups:
    tmp_genomes = reference_dist.index.intersection( source_groups.groups[isolation_source] )
    if len(tmp_genomes) < 5:
        continue

    print isolation_source
    sns.kdeplot( squareform( reference_dist.loc[tmp_genomes][tmp_genomes].values ), ax=ax, label=isolation_source )
ax.legend(loc=1)
fig.tight_layout()
fig.savefig('within_niche-evo_dists-reference.pdf')

for group in t3ss_groups:

    if not isfile('tree_puzzle/nucleotide/distances/%s.phy.dist.tab' %group):
        continue

    print group

    tmp_dist = pd.read_table('tree_puzzle/nucleotide/distances/%s.phy.dist.tab' %group, index_col=0)
    new_index = []
    for index in tmp_dist.index:
        new_index.append( index.split('|')[1] )

    tmp_dist.index   = new_index
    tmp_dist.columns = new_index

    fig, ax = plt.subplots()
    for isolation_source in source_groups.groups:
        tmp_genomes = tmp_dist.index.intersection( source_groups.groups[isolation_source] )
        if len(tmp_genomes) < 5:
            continue

        print '\t%s' %isolation_source
        sns.kdeplot( squareform( tmp_dist.loc[tmp_genomes][tmp_genomes].values ), label=isolation_source, ax=ax )
    ax.legend(loc=1)
    fig.tight_layout()
    fig.savefig('evo_dists_within_niche/within_niche-evo_dists-%s.pdf' %group)
    fig.clear()
    plt.close()
print 'yeah'

################################################################################

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
