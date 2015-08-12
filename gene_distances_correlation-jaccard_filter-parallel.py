from matplotlib import use
use('Agg')
from matplotlib import pyplot as plt
from functools import partial
from time import time
from pickle import load, dump
import MySQLdb as sql
import pandas as pd
from os import listdir, chdir, system
from scipy.stats import spearmanr, pearsonr
from scipy.spatial.distance import squareform, jaccard
from itertools import combinations, izip_longest
import multiprocessing
import networkx as nx
import numpy as np
from statsmodels.sandbox.stats.multicomp import fdrcorrection0 as fdr
from progressbar import ProgressBar

chdir('/Volumes/Macintosh HD 2/thiberio/virulence_factors/network-jaccard_filter')

groups = load( open( 'homologous_groups-merged.pkl' ) )
pb = pd.read_table('../../presence_absence-merged.tab', index_col=0, nrows=105)

######################################################
########################### finally, assess distances!
######################################################
nice_families = {}
for dist_table in listdir('distances'):

    if not dist_table.endswith('.tab'):
        continue

    group = dist_table.replace('.phy.dist.tab', '')
    group = group.replace('-', '&')

    if pb[group].sum() < 10:
        continue

    flag = False

    df = pd.read_table('distances/%s' %dist_table, index_col=0)

    if df.shape[0] == len(groups[group]):

        species = []
        for index in df.index:
            species.append( index.split('|')[1] )
        df.index   = species
        df.columns = species

        nice_families[group] = df.copy()
        continue
    else:
        pass

species_combinations = []
for pair in combinations(pb.index, 2):
    species_combinations.append( frozenset( pair ) )
all_species_distances = pd.DataFrame(index=nice_families.keys(), columns=species_combinations)

for group in nice_families.keys():

    df = nice_families[group].copy()

    indexes = []
    for pair in combinations(df.index, 2):
        indexes.append( frozenset(pair) )

    all_species_distances.loc[group][indexes] = squareform(df)
all_species_distances = all_species_distances.astype(float)

columns = all_species_distances.columns.tolist()
for n in range(len(columns)):
    columns[n] = '--'.join(columns[n])
all_species_distances.columns = columns

all_species_distances.to_csv('all_species_distances.tab', sep='\t')


######################################################
######################## finally, assess correlations!
######################################################

#
# if want to read a pre-computed distance DataFrame, uncomment line below
#all_species_distances = pd.read_table('all_species_distances.tab', index_col = 0)

ignored_species = ['A_schubertii_CECT4240T', 'A_diversa_CECT4254T', 'A_simiae_CIP107798T']
columns_to_remove = []
for column in all_species_distances.columns:
    species = column.split('--')
    if species[0] in ignored_species or species[1] in ignored_species:
        columns_to_remove.append(column)
all_species_distances.drop(columns_to_remove, axis=1, inplace=True)
pb.drop(ignored_species, inplace=True)

def assess_corr(df, pb, pairs):
    condensed_corr = []
    condensed_pval = []
    
    for pair in pairs:

        if jaccard(pb[pair[0]], pb[pair[1]]) > 0.2:
            (rho, p) = (np.nan, np.nan)
        else:
            tmp_df = df.loc[list(pair)].copy().dropna(axis=1, how='any').T
            (rho, p) = spearmanr(tmp_df)

        condensed_corr.append(rho)
        condensed_pval.append(p  )

    return (condensed_corr, condensed_pval)

group_pairs = list(combinations(all_species_distances.index, 2))
num_of_threads = 10
num_of_comparisons = len(group_pairs)
print "\t** Breaking datasets..."
avg = num_of_comparisons / float(num_of_threads)
datasets = []
last = 0.0
while last < len(group_pairs):
    datasets.append(group_pairs[int(last):int(last + avg)])
    last += avg
print "\t... done!\n"

print "\t** Getting shit done..."
start_time = time()
pool = multiprocessing.Pool(processes=num_of_threads)
func = partial(assess_corr, all_species_distances, pb)
results = pool.map(func, datasets)
pool.close()
pool.join()

print "\t... almost, just organizing ..."

condensed_corr = []
condensed_pval = []
for thread_result in results:
    condensed_corr.extend(thread_result[0])
    condensed_pval.extend(thread_result[1])
condensed_corr = np.array(condensed_corr)
condensed_pval = np.array(condensed_pval)
print "\t... done!\n"
print time() - start_time

out = open('condensed_corr.pkl', 'wb')
dump(condensed_corr, out)
out.close()
out = open('condensed_pval.pkl', 'wb')
dump(condensed_pval, out)
out.close()

#
# if wanna load pre-computed condensed correlation and p-value matrices, uncomment bellow
#print "\t**Loading condensed matrices..."
#condensed_corr = load( open( 'condensed_corr.pkl' ) )
#condensed_pval = load( open( 'condensed_pval.pkl' ) )
#print "\t... done!\n"

print "\t**Correcting p-values..."
pvals_tested = condensed_pval[pd.notnull(condensed_pval)]
(rejecteds, corrected_pval) = fdr(pvals_tested, alpha=0.01)

pos = 0
should_reject_rho = []
for uncorrect in condensed_corr:
    if np.isnan(uncorrect):
        should_reject_rho.append(False)
    else:
        should_reject_rho.append(rejecteds[pos])
        pos += 1
print "\t... done!\n"

print "\t**Organizing condensed matrices into something you can understand (you dumbass!)..."
uncorrected_corr_matrix = squareform(condensed_corr)
for n in range(uncorrected_corr_matrix.shape[0]):
    uncorrected_corr_matrix[n][n] = 1.
uncorrected_corr_df = pd.DataFrame(index=all_species_distances.index, columns=all_species_distances.index, data=uncorrected_corr_matrix)

rejecteds_df = pd.DataFrame(index=all_species_distances.index, columns=all_species_distances.index, data=squareform(should_reject_rho) > 0)
corr_df = uncorrected_corr_df[rejecteds_df]
corr_df.to_csv('significant_gene_families_dist_correlations.tab', sep='\t')
print "\t... done!\n"

#
# network time!
#
graph = nx.Graph()

indexes = corr_df.index
graph.add_nodes_from(corr_df.index)

print '\t**Adding edges ...'
edges_between = list( combinations( indexes, 2 ) )
weights       = squareform( corr_df.fillna(0).values )
for n in range( len( weights) ):
    edge_weight = weights[n]
    if edge_weight >= 0.9:
        connected_nodes = edges_between[n]
        graph.add_edge( connected_nodes[0], connected_nodes[1], weight = edge_weight )

degrees = graph.degree()
nodes_with_no_connections = []
for node in degrees.iteritems():
    if node[1] == 0:
        nodes_with_no_connections.append(node[0])
graph.remove_nodes_from(nodes_with_no_connections)

group_length = []
node_sizes   = []
node_colors  = []
for node in graph.nodes():
    node_sizes.append( len( groups[node] ) )
    node_colors.append( degrees[node] )

    phylip_header = open('phylip/%s.phy' %node.replace('&', '-')).xreadlines().next()
    group_length.append( int( phylip_header.split()[1] ) )
print '... done!'

print '\t**Generating spring layout ...'
graph_layout = nx.spring_layout(graph, iterations=500)
print '... done!'

plt.figure()
plt.xticks([])
plt.yticks([])

print '\t**Drawing nodes ...'
nx.draw_networkx_nodes(graph, pos=graph_layout, node_size=node_sizes, node_color=node_colors, with_labels=False, alpha=0.75, cmap=plt.get_cmap('Blues'))
print '\t**Saving pdf (no edges) ...'
plt.savefig('correlation_network-ignored_species-no_edges.pdf', bbox_inches='tight')#, dpi=100, figsize=(15,15))
print '\t**Drawing edges ...'
nx.draw_networkx_edges(graph, pos=graph_layout, with_labels=False, alpha=0.3)
print '\t**Saving final pdf ...'
plt.savefig('correlation_network-ignored_species.pdf', bbox_inches='tight')#, dpi=100, figsize=(15,15))
print '... done!'
plt.close()
######################################################
############### test if gene length influences degree!
######################################################
plt.figure()
plt.hist(group_length, bins=100, alpha=0.7, edgecolor='white')
plt.tight_layout()
plt.savefig('hist_length.pdf')
plt.close()

plt.figure()
plt.hist(node_colors, bins=100, alpha=0.7, edgecolor='white')
plt.tight_layout()
plt.savefig('hist_degree.pdf')
plt.close()

plt.figure()
plt.hist(node_sizes, bins=100, alpha=0.7, edgecolor='white')
plt.tight_layout()
plt.savefig('hist_num_of_species.pdf')
plt.close()

plt.figure()
plt.scatter(group_length, node_colors, alpha=0.3, c='g', edgecolor='none')
plt.xlim(xmin=0)
plt.ylim(ymin=0)
plt.xlabel('Gene family sequence length (MSA length)')
plt.ylabel('Gene family degree')
plt.tight_layout()
plt.savefig('length_VS_degree.pdf')
plt.close()

plt.figure()
plt.scatter(node_sizes, node_colors, alpha=0.3, c='g', edgecolor='none')
plt.xlabel('Num os species where gene family is present')
plt.ylabel('Gene family degree')
plt.xlim(xmin=0)
plt.ylim(ymin=0)
plt.tight_layout()
plt.savefig('num_of_sp_VS_degree.pdf')
plt.close()
######################################################

######################################################
############################################ MCL time!
######################################################
print '\t**Saving network matrix ...'
nx.write_weighted_edgelist(graph, 'mcl_input.abc')

print '\t**Writing in MCL-happy format ...'
system('mcxload -abc mcl_input.abc --stream-mirror -o mcl_input.mci -write-tab mcl_input.tab')

print '\t**Running MCL with different inflation values (-I) ...'
system('cat inflation_values | xargs -n1 -P10 mcl mcl_input.mci -I ')
system('clm dist --chain out.mcl_input.mci.I{14,20,25,30,35,40,45,50,55,60,65,70} > clm_dist')

print '\t**Running MCL with selected inflation ...'
system('mcl mcl_input.mci -I 2 -te 10 -o mcl_output')

system('mcxdump -icl mcl_output -tabr mcl_input.tab -o gene_clusters')
print '... done!'

gene_clusters = []
for line in open('gene_clusters').xreadlines():
    gene_clusters.append(line.split())
    print len(line.split())
print len(gene_clusters)

hm21 = 'A_veronii_Hm21'
hm21_mls = {
'dnaK' : '99.2089',
'dnaJ' : '99.2088',
'gltA' : '99.1610',
'groL' : '99.2695',
'gyrA' : '99.1403',
'gyrB' : '99.3369',
'recA' : '99.2915',
'rpoB' : '99.3720',
'rpoD' : '99.2705',
'dnaX' : '99.1560',
'radA' : '99.2751',
'zipA' : '99.2310',
'tsf'  : '99.2388',
'atpD' : '99.3397'
}
mls_groups = {}
for group in groups:
    if hm21 not in groups[group]:
        continue
    for gene in hm21_mls:
        if hm21_mls[gene] in groups[group][hm21]:
            mls_groups[gene] = group
            break

for cluster in gene_clusters:
    for group in mls_groups.values():
        if group in cluster:
            print group
            print len(cluster)
            print ''

for genes in combinations(mls_groups.keys(), 2):
    print genes

    pair = [mls_groups[genes[0]], mls_groups[genes[1]]]
    tmp = distances_df.loc[pair].T.copy()
    tmp.dropna(axis=0, how='any')

    (rho, pval) = spearmanr(tmp)

    plt.figure()
    tmp.plot(x=0, y=1, kind='scatter', alpha=0.5, edgecolor='none', c='k')
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)
    plt.xlabel(genes[0])
    plt.ylabel(genes[1])
    plt.title('rho: %.2f | p-value: %.2e' %(rho, pval))
    plt.savefig('mls_gene_corr/%s.pdf' %'-'.join(genes), bbox_inches='tight', dpi=300, figsize=(15,15))
    plt.close()
