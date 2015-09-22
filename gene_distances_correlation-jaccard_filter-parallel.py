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
from skbio.stats.distance import mantel, DistanceMatrix as dm

main_dir = '/Volumes/Macintosh HD 2/thiberio/virulence_factors/'

groups = load( open( '%s/homologous_groups-merged.pkl' %main_dir ) )
pb = pd.read_table('%s/presence_absence-merged.tab' %main_dir, index_col=0, nrows=105)
print 'yeah'

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
ignored_species = ['A_schubertii_CECT4240T', 'A_diversa_CECT4254T', 'A_simiae_CIP107798T']

dist_matrices = []
group_ids     = []
for n in listdir('distances'):
    if not n.endswith('.phy.dist.tab'):
        continue

    group_name = n.replace('.phy.dist.tab', '')
    tmp = pd.read_table('distances/%s' %n, index_col=0)
    taxa = []
    for taxon in tmp.index:
        taxa.append(taxon.split('|')[1])

    if len(taxa) > len( set( taxa ) ):
        continue

    tmp.index   = taxa
    tmp.columns = taxa

    should_be_dropped = list(np.intersect1d(taxa, ignored_species))

    if should_be_dropped:
        tmp.drop(should_be_dropped, axis=0, inplace=True)
        tmp.drop(should_be_dropped, axis=1, inplace=True)

    dist_matrices.append(  dm(tmp.as_matrix(), tmp.index) )
    group_ids.append(                          group_name )


def assess_corr(pairs):
    condensed_corr = []
    condensed_pval = []
    
    for pair in pairs:

        j_index =  len( np.intersect1d( pair[0].ids, pair[1].ids ) ) / float( len( np.union1d( pair[0].ids, pair[1].ids ) ) ) 
        if j_index < 0.8:
            (rho, p) = (np.nan, np.nan)
        else:
            (rho, p, n) = mantel(pair[0], pair[1], strict=False)

        condensed_corr.append(rho)
        condensed_pval.append(p  )

    return (condensed_corr, condensed_pval)

group_pairs = list(combinations(dist_matrices, 2))
num_of_threads = 15
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
results = pool.map(assess_corr, datasets)
pool.close()
pool.join()
print "\t... done!\n"
print time() - start_time

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





#
# if want to read a pre-computed distance DataFrame, uncomment line below
#all_species_distances = pd.read_table('all_species_distances.tab', index_col = 0)

ignored_species = ['A_schubertii_CECT4240T', 'A_diversa_CECT4254T', 'A_simiae_CIP107798T']
columns_to_remove = []
for column in all_species_distances.columns:
    species = column.split('--')
    if species[0] in ignored_species or species[1] in 
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
uncorrected_corr_df = pd.DataFrame(index=group_ids, columns=group_ids, data=uncorrected_corr_matrix)

rejecteds_df = pd.DataFrame(index=group_ids, columns=group_ids, data=squareform(should_reject_rho) > 0)
corr_df = uncorrected_corr_df[rejecteds_df]
corr_df.to_csv('significant_gene_families_mantel_dist_correlations.tab', sep='\t')
print "\t... done!\n"

#
# network time!
#
corr_df = pd.read_table('significant_gene_families_mantel_dist_correlations.tab', index_col=0)
graph = nx.Graph()

indexes = corr_df.index
graph.add_nodes_from(corr_df.index)

print '\t**Adding edges ...'
edges_between = list( combinations( indexes, 2 ) )
weights       = squareform( corr_df.fillna(0).values )
threshold = 0.7
for n in range( len( weights) ):
    edge_weight = weights[n]
    if edge_weight >= threshold:
        connected_nodes = edges_between[n]
        graph.add_edge( connected_nodes[0], connected_nodes[1], weight = edge_weight )

degrees = graph.degree()
nodes_with_no_connections = []
for node in degrees.iteritems():
    if node[1] == 0 or node[0] not in groups:
        nodes_with_no_connections.append(node[0])
graph.remove_nodes_from(nodes_with_no_connections)

group_length = []
node_sizes   = []
node_colors  = []
for node in graph.nodes():
    node_sizes.append( len( groups[node] ) )
    node_colors.append( degrees[node] )

#    phylip_header = open('phylip/%s.phy' %node.replace('&', '-')).xreadlines().next()
#    group_length.append( int( phylip_header.split()[1] ) )
print '... done!'

print '\t**Generating spring layout ...'
graph_layout = nx.spring_layout(graph)#, iterations=500)
print '... done!'

plt.figure()
plt.xticks([])
plt.yticks([])

print '\t**Drawing nodes ...'
nx.draw_networkx_nodes(graph, pos=graph_layout, node_size=node_sizes, node_color=node_colors, with_labels=False, alpha=0.75, cmap=plt.get_cmap('Blues'))
print '\t**Saving pdf (no edges) ...'
plt.savefig('mantel_correlation_network-no_cutoff-ignored_species-no_edges.pdf', bbox_inches='tight')#, dpi=100, figsize=(15,15))
print '\t**Drawing edges ...'
nx.draw_networkx_edges(graph, edgelist=edgelist, pos=graph_layout, with_labels=False, alpha=0.3)
print '\t**Saving final pdf ...'
plt.savefig('mantel_correlation_network-no_cutoff-ignored_species.pdf', bbox_inches='tight')#, dpi=100, figsize=(15,15))
print '... done!'
plt.close()

######################################################
############################################ MCL time!
######################################################
print '\t**Saving network matrix ...'
nx.write_weighted_edgelist(graph, 'mcl_mantel_input-%.1f.abc' %threshold)

print '\t**Writing in MCL-happy format ...'
system('mcxload -abc mcl_mantel_input-%.1f.abc --stream-mirror -o mcl_mantel_input-%.1f.mci -write-tab mcl_mantel_input-%.1f.tab' %(threshold, threshold, threshold))

print '\t**Running MCL with different inflation values (-I) ...'
system('cat inflation_values | xargs -n1 -P10 mcl mcl_mantel_input-%.1f.mci -I ' %threshold)
system('clm dist --chain out.mcl_mantel_input-%.1f.mci.I{14,20,25,30,35,40,45,50,55,60,65,70} > clm_dist-%.1f' %threshold)

print '\t**Running MCL with selected inflation ...'
system('mcl mcl_mantel_input-%.1f.mci -I 4.5 -te 10 -o mcl_mantel_output-%.1f' %(threshold, threshold))

system('mcxdump -icl mcl_mantel_output-%.1f -tabr mcl_mantel_input-%.1f.tab -o gene_clusters-I45-mantel-%.1f' %(threshold, threshold, threshold))
print '... done!'

gene_clusters = []
for line in open('test/gene_clusters-0.7').xreadlines():
    gene_clusters.append(line.split())
top20 = []
for cluster in gene_clusters[:20]:
    top20.extend(cluster)
top40 = []
for cluster in gene_clusters[:40]:
    top40.extend(cluster)

weights = {}
nodes = []
for cluster in range(10):
    nodes.extend(gene_clusters[cluster])
for edge in graph.subgraph(nodes).edges(data=True):
    (a,b) = (False, False)
    for cluster in range(10):
        if edge[0] in gene_clusters[cluster]:
            a = cluster
        if edge[1] in gene_clusters[cluster]:
            b = cluster
        if a and b:
            break

    keys = [a,b]
    keys.sort()
    key = str(keys[0])+','+str(keys[1])
    if key not in weights:
        weights[key] = []

    weights[key].append(edge[2]['weight'])

print 'yeah'

for i in range(10):
    i = str(i)
    
    for j in range(10):
        j = str(j)

        if i == j:
            continue

        if '%s,%s' %(i,j) in weights:
            t,p = ttest_ind(weights['%s,%s' %(i,i)], weights['%s,%s' %(i,j)])
            print '%s vs %s-> t: %.2f\tp:%.2e' %(i, j, t, p)

tree_comp_dataset = {}
for i in range(10):

    for j in range(10):

        pairs = []
    
        if i < j:
            break
        elif i == j:
            pairs = list( combinations( gene_clusters[i],2 ) )
        else:
            pairs = list( product( gene_clusters[i], gene_clusters[j] ) )

        print '%i vs %i' %(i,j)

        if not pairs:
            continue

        shuffle(pairs)

        key = '%i,%i' %(i, j)
        if key not in tree_comp_dataset:
            tree_comp_dataset[key] = []

        for pair in pairs:
            j_dist =  jaccard(pb[pair[0]], pb[pair[1]])

            if j_dist > 0.2:
                continue

            tree_comp_dataset[key].append(pair)
            if len(tree_comp_dataset[key]) == 100:
                break

        if len(tree_comp_dataset[key]) < 20:
            print '\t%s is empty! =/' %key
            tree_comp_dataset.pop(key)
print 'yeah'

def run_clann_rfdists(cluster_pair, tree_pairs):
    tmp_input_file_name = 'tmp_input_rf_%f' %random.random()
    tmp_output_file_name = 'tmp_output_rf_%f' %random.random()

    for pair in tree_pairs:
        system('cat fasttree-no_gene_id/%s.tre fasttree-no_gene_id/%s.tre > %s' %(pair[0], pair[1], tmp_input_file_name))
        system('echo "execute %s; rfdists filename=%s; quit;" | clann' %(tmp_input_file_name, tmp_output_file_name))
        system('cat %s >> %s.rf' %(tmp_output_file_name, cluster_pair))

    system('rm %s' %tmp_input_file_name)
    system('rm %s' %tmp_output_file_name)
map(lambda x: run_clann_rfdists(x, tree_comp_dataset[x]), tree_comp_dataset.keys())

def run_dp_wrfdist(tree_pairs):
    dists = []

    for pair in tree_pairs:
        t1 = tree.get_from_path('fasttree-no_gene_id/%s.tre' %pair[0], 'newick', taxon_namespace=tns)
        t2 = tree.get_from_path('fasttree-no_gene_id/%s.tre' %pair[1], 'newick', taxon_namespace=tns)
        
        t1.encode_bipartitions()
        t2.encode_bipartitions()

        dists.append( treecompare.weighted_robinson_foulds_distance( t1, t2 ) )

    return dists
yeah = {}
yeah = map(lambda x: run_dp_wrfdist(tree_comp_dataset[x]), tree_comp_dataset.keys())

rfdists = {}
for cluster_pair in tree_comp_dataset:
    rfdists[cluster_pair] = []
    for dist in open('%s.rf' %cluster_pair).read().split():
        rfdists[cluster_pair].append(float(dist))

for pair in rfdists:
    if pair[0] == pair[2]:
        continue

    t,p = ttest_ind(rfdists[pair], rfdists['%s,%s' %(pair[0],pair[0])])
    print '%s vs %s,%s-> t: %.2f\tp:%.2e' %(pair, pair[0], pair[0], t, p)
    t,p = ttest_ind(rfdists[pair], rfdists['%s,%s' %(pair[2],pair[2])])
    print '%s vs %s,%s-> t: %.2f\tp:%.2e' %(pair, pair[2], pair[2], t, p)
print 'yeah'




        

for cluster in gene_clusters:
    for group in mls_groups:
        if mls_groups[group] in cluster:
            print group
            print gene_clusters.index(cluster)
            print ''

#
# network time, again!
#
graph = nx.Graph()

indexes = corr_df.index
graph.add_nodes_from(corr_df.index)

print '\t**Adding edges ...'
edges_between = list( combinations( indexes, 2 ) )
weights       = squareform( corr_df.fillna(0).values )
for n in range( len( weights) ):
    edge_weight = weights[n]
    if edge_weight >= threshold:
        connected_nodes = edges_between[n]
        graph.add_edge( connected_nodes[0], connected_nodes[1], weight = edge_weight )

degrees = graph.degree()
nodes_with_no_connections = []
for node in degrees.iteritems():
    if node[1] == 0 or node[0] not in groups or node[0] not in top20:
        nodes_with_no_connections.append(node[0])
graph.remove_nodes_from(nodes_with_no_connections)

node_sizes   = []
node_colors  = []
for node in graph.nodes():
    node_sizes.append( len( groups[node] ) )
    flag = True
    for pos in range(20):
        if node in gene_clusters[pos]:
            node_colors.append(pos)
            flag = False
            break
    if flag:
        continue
        node_colors.append(0)
print '... done!'

print '\t**Generating spring layout ...'
graph_layout = nx.spring_layout(graph)
print '... done!'

cNorm  = colors.Normalize(vmin=0, vmax=41)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plt.get_cmap('Set1'))

plt.figure()
plt.xticks([])
plt.yticks([])

print '\t**Drawing nodes ...'
nx.draw_networkx_nodes(graph, nodelist=not_top40, pos=graph_layout, node_size=node_sizes, node_color=['w']*len(not_top40), label='not top 40', with_labels=False, alpha=0.3)
for cluster in range(40):
    nx.draw_networkx_nodes(graph, nodelist= gene_clusters[cluster], pos=graph_layout, node_size=node_sizes, node_color=[scalarMap.to_rgba(cluster)] * len(gene_clusters[cluster]), label='cluster_%i' %cluster, with_labels=False, alpha=0.6)
plt.legend(loc=9, scatterpoints=1, bbox_to_anchor=(0.5,0), ncol=9, fontsize='xx-small', columnspacing=0.3, handletextpad=0)
print '\t**Saving pdf (no edges) ...'
plt.savefig('xmantel_correlation_network-uhuhuhtop40-I45-%.1f-ignored_species-no_edges10.pdf' %threshold, bbox_inches='tight')#, dpi=100, figsize=(15,15))
print '\t**Drawing edges ...'
nx.draw_networkx_edges(graph, edgelist=edgelist, pos=graph_layout, with_labels=False, alpha=0.3)
print '\t**Saving final pdf ...'
plt.savefig('xmantel_correlation_network-%.1f-ignored_species1.pdf' %threshold, bbox_inches='tight')#, dpi=100, figsize=(15,15))
print '... done!'
plt.close()

for cluster in gene_clusters[:20]:
    print gene_clusters.index(cluster)
    print nx.average_clustering(graph.subgraph(cluster))
    print ''







for cluster in gene_clusters:
    if len(cluster) < 10:
        break
    system('cat fasttree-no_gene_id/%s.tre > gene_cluster-%i.tre' %( '.tre fasttree-no_gene_id/'.join(cluster) , gene_clusters.index( cluster ) ) )
    

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
'atpD' : '99.3397',
'metG' : '99.1415'
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
    for group in mls_groups:
        if mls_groups[group] in cluster:
            print group
            print gene_clusters.index(cluster)
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



def strip_gene_name(tree_name):
    tree = read('fasttree/%s' %tree_name, 'newick')
    for leaf in tree.get_terminals():
        leaf.name = leaf.name.split('|')[1]
    write(tree, 'fasttree-no_gene_id/%s' %tree_name, 'newick')
from multiprocessing import Pool
pool = Pool(processes=10)
pool.map(strip_gene_name, listdir('fasttree/'))
pool.close()
pool.join()
