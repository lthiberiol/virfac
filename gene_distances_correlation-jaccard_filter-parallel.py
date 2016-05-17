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
from scipy import stats
from scipy.spatial.distance import squareform, jaccard
from itertools import combinations, izip_longest
import multiprocessing
import networkx as nx
import numpy as np
from statsmodels.sandbox.stats.multicomp import fdrcorrection0 as fdr
from skbio.stats.distance import mantel, DistanceMatrix as dm
from random import sample

main_dir = '/Volumes/Macintosh HD 2/thiberio/virulence_factors/'

groups = load( open( '%s/homologous_groups-merged.pkl' %main_dir ) )
pb = pd.read_table('%s/presence_absence-merged.tab' %main_dir, index_col=0, nrows=105)
print 'yeah'

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

    if len( taxa ) != len( set( taxa ) ):
        continue

    tmp.index   = taxa
    tmp.columns = taxa

    should_be_dropped = list(np.intersect1d(taxa, ignored_species))

    if should_be_dropped:
        tmp.drop(should_be_dropped, axis=0, inplace=True)
        tmp.drop(should_be_dropped, axis=1, inplace=True)

    dist_matrices.append(  dm(tmp.as_matrix(), tmp.index) )
    group_ids.append(                          group_name )

def is_intercept_on_origin(pairs, alpha=0.05):
    intercepts = []

    for dm1, dm2 in pairs:

        intersection = np.intersect1d( dm1.ids, dm2.ids )
        union        = np.union1d(     dm1.ids, dm2.ids )
        j_index =  len( intersection ) / float( len( union ) ) 

        if j_index < 0.3:
            intercepts.append(None)
            continue

        x = dm1.filter(intersection).condensed_form()
        y = dm2.filter(intersection).condensed_form()

        n = len(x)             # number of samples     
        Sxx = np.sum(x**2) - np.sum(x)**2/n
        Syy = np.sum(y**2) - np.sum(y)**2/n    # not needed here
        Sxy = np.sum(x*y) - np.sum(x)*np.sum(y)/n    

        (slope, intercept, rcoef, pvalue, slope_std_err) = stats.linregress(x, y)

        # Residuals
        fit = lambda xx: intercept + slope*xx    
        residuals = y - fit(x)
        var_res = np.sum(residuals**2)/(n-2)
        sd_res = np.sqrt(var_res)

        # Confidence intervals
        se_slope = sd_res/np.sqrt(Sxx)
        se_intercept = sd_res*np.sqrt(np.sum(x**2)/(n*Sxx))

        df = n-2                            # degrees of freedom
        tval = stats.t.isf(0.05/2., df)     # appropriate t value

        ci_intercept = intercept + tval*se_intercept*np.array([-1,1])

        intercepts.append( True if ci_intercept[0] <= 0 <= ci_intercept[1] else False )

    return intercepts

def assess_corr(pairs):
    condensed_corr = []
    condensed_pval = []
    
    for pair in pairs:

        intersection = np.intersect1d( pair[0].ids, pair[1].ids )
        union        = np.union1d(     pair[0].ids, pair[1].ids )
        j_index =  len( intersection ) / float( len( union ) ) 

        if j_index < 0.8:
            (rho, p) = (np.nan, np.nan)
        else:
            (rho, p, n) = mantel(pair[0], pair[1], strict=False)
            intercept_flag = is_intercept_on_origin( pair[0].filter(intersection), pair[1].filter(intersection) )

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
results = pool.map(is_intercept_on_origin, datasets)
pool.close()
pool.join()
print "\t... done!\n"
print time() - start_time

print "\t... almost, just organizing ..."

condensed_intercept_on_origin = []
for thread_result in results:
    condensed_intercept_on_origin.extend(thread_result)
condensed_intercept_on_origin = np.array(condensed_intercept_on_origin)
print "\t... done!\n"
intercept_df = pd.DataFrame(index=group_ids, columns=group_ids, data=squareform(condensed_intercept_on_origin))

condensed_corr = []
condensed_pval = []
for thread_result in results:
    condensed_corr.extend(thread_result[0])
    condensed_pval.extend(thread_result[1])
condensed_corr = np.array(condensed_corr)
condensed_pval = np.array(condensed_pval)
print "\t... done!\n"
print time() - start_time

out = open('condensed_mantel_corr.pkl', 'wb')
dump(condensed_corr, out)
out.close()
out = open('condensed_mantel_uncorrected_pval.pkl', 'wb')
dump(condensed_pval, out)
out.close()
print "\t... done!\n"

#
# if wanna load pre-computed condensed correlation and p-value matrices, uncomment bellow
#print "\t**Loading condensed matrices..."
#condensed_corr = load( open( 'condensed_mantel_corr.pkl' ) )
#condensed_pval = load( open( 'condensed_mantel_uncorrected_pval.pkl' ) )
#print "\t... done!\n"

print "\t**Correcting p-values..."
pvals_tested = condensed_pval[pd.notnull(condensed_pval)]
(rejecteds, corrected_pval) = fdr(pvals_tested, alpha=0.05)

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
xcorr_df = uncorrected_corr_df[rejecteds_df]
corr_df.to_csv('significant_gene_families_mantel_dist_correlations.tab', sep='\t')
print "\t... done!\n"


fm_df = pd.read_table('fitch-margoliash_distances.tab', index_col=0)
condensed_fm = squareform(fm_df.fillna(10))

hell1 = []
yeah1 = []
for x,y,z in zip(condensed_fm, condensed_intercept_on_origin, squareform(corr_df.fillna(0))):
    if x<10 and z:
        hell1.append(x)
        yeah1.append(z)

contig1_distances = pd.read_table('/Users/Thiberio/work/A_veronii_Hm21-conversion/hm21_contig1_genomic_distances.tab', index_col=0)
contig1_distances.index = contig1_distances.columns
contig2_distances = pd.read_table('/Users/Thiberio/work/A_veronii_Hm21-conversion/hm21_contig2_genomic_distances.tab', index_col=0)
contig2_distances.index = contig2_distances.columns
adjacent_pairs = load(open('hm21_adjacent_pairs.pkl'))
cect839_genes = set()
for pair in adjacent_pairs:
    cect839_genes.update(pair)

gene2group = {}
group2gene = {}
#species = 'A_hydrophila_CECT839T'
species = 'A_veronii_Hm21'
for group in groups:
    if species not in groups[group]:
        continue
    if group not in fm_df.index:
        continue
    gene = groups[group][species][0]
    if gene not in cect839_genes:
        continue
    group2gene[group] = gene
    gene2group[gene] = group

clusters = []
for line in open('gene_clusters-0.7-I55').xreadlines():
    clusters.append(set(line.split()))

adjacent_pairs = load(open('hm21_adjacent_pairs.pkl'))
adjacent_groups = set()
for pair in adjacent_pairs:
    pair = list(pair)
    if pair[1] in gene2group and pair[0] in gene2group:
        adjacent_groups.add( frozenset([gene2group[pair[0]], gene2group[pair[1]]]) )

xablau = clusters[6].intersection(group2gene.keys())
observed_pairs = list( combinations(xablau,2) )
for n in range(len(observed_pairs)):
    observed_pairs[n] = frozenset(observed_pairs[n])
len( adjacent_groups.intersection(observed_pairs) )


manager = multiprocessing.Manager()
yeah = manager.list(adjacent_groups)
def random_tries(random_sample):
    random_test = set()
    for n in combinations(random_sample, 2):
        random_test.add( frozenset(n) )
    return len( random_test.intersection(adjacent_groups) )

random_samples = []
for n in range(1000):
    random_samples.append( sample(group2gene.keys(), len(xablau)))
    
pool = multiprocessing.Pool(processes=10)
random_intersection = pool.map(random_tries, random_samples)
pool.close()
pool.join()
print 'yeah'

ser = pd.Series(random_intersection)
plt.figure()
ser.plot(kind='density')
plt.plot(len( adjacent_groups.intersection(observed_pairs) ), 0.001, 'ro', alpha=0.6, markersize=10)
plt.savefig('random_adjacent_pairs-cect839-cluster6.pdf', dpi=300, figsize=(15, 15), bbox_inches='tight')

out = open('random_adjacent_pairs-cluster1.pkl', 'wb')
dump(random_intersection, out)
out.close()

hm21_pathways = load(open('hm21_pathways.pkl'))
pathways_fixed_order = set()
for pathways in hm21_pathways.values():
    pathways_fixed_order.update(pathways)
pathways_fixed_order = tuple(pathways_fixed_order)

reference_count = pd.Series(data=0.0, index=pathways_fixed_order)
xablau = clusters[29].intersection(group2gene.keys())
for group in xablau:
    if group not in group2gene:
        continue
    
    gene = group2gene[group]
    
    if gene not in hm21_pathways:
        continue

    for pathway in hm21_pathways[gene]:
        reference_count[pathway] += 1

random_samples = pd.DataFrame(columns = pathways_fixed_order)
for n in range(1000):
    random_sample= sample(group2gene.keys(), len(xablau))
    tmp_count = pd.Series(data=0.0, index=pathways_fixed_order)

    for group in random_sample:
        gene = group2gene[group]
        
        if gene not in hm21_pathways:
            continue

        for pathway in hm21_pathways[gene]:
            tmp_count[pathway] += 1
    random_samples.loc[n] = tmp_count
print 'yeah'

non_randomly_enriched = reference_count > random_samples.mean() + random_samples.std()*1
non_randomly_enriched = reference_count[non_randomly_enriched].index

ind = np.arange(len(reference_count[non_randomly_enriched]))
width = 0.35
xticks = []
for n in non_randomly_enriched:
    xticks.append(n.replace('path:', ''))

plt.figure()
plt.bar(ind, reference_count[non_randomly_enriched], width, color='g', edgecolor='none', alpha=0.6, label='Pathway count')
plt.bar(ind+width, random_samples[non_randomly_enriched].mean(), width, color='b', edgecolor='none', alpha=0.6, yerr=random_samples[non_randomly_enriched].std(), label='Random counts')
plt.xticks(ind+0.5, xticks, rotation=45, ha='right')
plt.ylim(ymin=0)
plt.legend()
plt.savefig('pathway_count-hm21-cluster29-rand1000.pdf', bbox_inches='tight', figsize=(300,10), dpi=300)
plt.close()
print 'yeah'


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
threshold = 0.95
for n in range( len( weights) ):
    edge_weight = weights[n]
    if edge_weight > threshold:
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

edgelist = []
neighbor_count = {}
for edge in graph.edges(data=True):
    if edge[2]['weight'] >= 0.9:
        edgelist.append(edge)
        for node in edge[:2]:
            if node not in neighbor_count:
                neighbor_count[node] = 0.0
            neighbor_count[node] += 1.0

sorted_nodes_by_degree = sorted(neighbor_count, key=neighbor_count.__getitem__)
sorted_nodes_by_degree.reverse()

#    phylip_header = open('phylip/%s.phy' %node.replace('&', '-')).xreadlines().next()
#    group_length.append( int( phylip_header.split()[1] ) )
print '... done!'

out = open('cluster0_nodes_sorted_by_degree', 'wb')
for node in sorted_nodes_by_degree:
    out.write( '%s: %i\n' %( node, neighbor_count[node] ) )
out.close()


print '\t**Generating spring layout ...'
graph_layout = nx.spring_layout(graph)#, iterations=500)
print '... done!'

plt.figure()
plt.xticks([])
plt.yticks([])

print '\t**Drawing nodes ...'
nx.draw_networkx_nodes(graph, pos=graph_layout, node_size=node_sizes, node_color=node_colors, with_labels=False, alpha=0.75, cmap=plt.get_cmap('Blues'))
print '\t**Saving pdf (no edges) ...'
plt.savefig('mantel_correlation_network-cluster0-no_edges.pdf', bbox_inches='tight', dpi=100, figsize=(15,15))
print '\t**Drawing edges ...'
nx.draw_networkx_edges(graph, edgelist=edgelist, pos=graph_layout, with_labels=False, alpha=0.3)
print '\t**Saving final pdf ...'
plt.savefig('mantel_correlation_network-cluster0.pdf', bbox_inches='tight')#, dpi=100, figsize=(15,15))
print '... done!'
plt.close()

######################################################
#################### Compare correlations with rfdists
######################################################
def run_clann_rfdists(tree_pairs):
    distances = []
    valid_pairs=[]

    tmp_input_file_name = 'tmp_input_rf_%f' %random.random()
    tmp_output_file_name = 'tmp_output_rf_%f' %random.random()

    for pair in tree_pairs:
        system('cat fasttree-no_gene_id/%s.tre fasttree-no_gene_id/%s.tre > %s' %(pair[0], pair[1], tmp_input_file_name))
        if not system('echo "execute %s; rfdists filename=%s; quit;" | clann' %(tmp_input_file_name, tmp_output_file_name)):
            distances.append( float( open( tmp_output_file_name ).read().strip() ) )
            valid_pairs.append(True)
        else:
            valid_pairs.append(False)

#    system('rm %s' %tmp_input_file_name)
#    system('rm %s' %tmp_output_file_name)

    return (valid_pairs, distances)

all_pairs = np.array( list( combinations( corr_df.index, 2 ) ) )
all_corrs = squareform( corr_df.fillna(0) )

condensed_corrs  = all_corrs[all_corrs != 0]
correlated_pairs = all_pairs[all_corrs != 0]

test_positions = sample(range( len( condensed_corrs ) ), 1000) 
test_pairs = map( lambda x: correlated_pairs[x], test_positions)
print 'yeah'

num_of_threads = 10
num_of_comparisons = len(test_pairs)
avg = num_of_comparisons / float(num_of_threads)
datasets = []
last = 0.0
while last < num_of_comparisons:
    datasets.append(test_pairs[int(last):int(last + avg)])
    last += avg


pool = Pool(processes=10)
results = pool.map(run_clann_rfdists, datasets)
pool.close()
pool.join()

valid_pairs = []
tree_dists  = []
for thread in results:
    valid_pairs.extend(thread[0])
    tree_dists.extend( thread[1])

equivalent_corrs = []
yeah = np.array(test_positions)
for pos in yeah[valid_pairs]:
    equivalent_corrs.append(condensed_corrs[pos])

plt.figure()
plt.scatter(equivalent_corrs, tree_dist)
plt.savefig('corr_vs_rfdist.pdf')

######################################################
############################################ MCL time!
######################################################
print '\t**Saving network matrix ...'
nx.write_weighted_edgelist(graph, 'mantel_network.abc')

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
for line in open('gene_clusters-p0.05-0.7-I25').xreadlines():
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

def run_95_rfdists(cluster_pair):
    for pair in b[cluster_pair]:
        system('python robinson_foulds-pw.py fasttree-no_gene_id/%s.tre fasttree-no_gene_id/%s.tre >> %s.rf' %(pair[0], pair[1], cluster_pair))
pool = multiprocessing.Pool(processes=15)
manager = multiprocessing.Manager()
b = manager.dict( deepcopy( tree_comp_dataset ) )
pool.map(run_95_rfdists, b.keys())

rfdists95 = {}
for cluster_pair in tree_comp_dataset:
    rfdists95[cluster_pair] = []
    for dist in re.findall('^\(\S+,\s(\S+)\)$' ,open('%s.rf' %cluster_pair).read(), re.M):
        rfdists95[cluster_pair].append(float(dist))

for pair in rfdists95:
    if pair[0] == pair[2]:
        continue

    t,p = ttest_ind(rfdists95[pair], rfdists95['%s,%s' %(pair[0],pair[0])])
    print '%s vs %s,%s-> t: %.2f\tp:%.2e' %(pair, pair[0], pair[0], t, p)
    t,p = ttest_ind(rfdists95[pair], rfdists95['%s,%s' %(pair[2],pair[2])])
    print '%s vs %s,%s-> t: %.2f\tp:%.2e' %(pair, pair[2], pair[2], t, p)

def run_clann_rfdists(tree_pairs):
    distances = []
    tmp_input_file_name = 'tmp_input_rf_%f' %random.random()
    tmp_output_file_name = 'tmp_output_rf_%f' %random.random()

    for pair in tree_pairs:
        system('cat fasttree-no_gene_id/%s.tre fasttree-no_gene_id/%s.tre > %s' %(pair[0], pair[1], tmp_input_file_name))
        system('echo "execute %s; rfdists filename=%s; quit;" | clann' %(tmp_input_file_name, tmp_output_file_name))
        distances.append( float( read( 'cat %s' %tmp_output_file_name ).open().strip() ) )

    return distances

correlated_pairs = list( combinations( corr_df.index, 2 ) )
correlations = squareform(corr_df)

map(lambda x: run_clann_rfdists(x, tree_comp_dataset[x]), tree_comp_dataset.keys())

rfdists = {}
for cluster_pair in tree_comp_dataset:
    rfdists[cluster_pair] = []
    for dist in open('%s-clann.rf' %cluster_pair).read().split():
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
