from Bio import Phylo as phy
import pandas as pd
from os import listdir, chdir
from itertools import combinations
from copy import deepcopy
import re
from multiprocessing import Pool
from random import sample
from time import time
from pickle import load, dump
from matplotlib import use
use('Agg')
from matplotlib import pyplot as plt

trees = phy.parse('RAxML_bootstrap.yeah', 'newick')
genomes = map( lambda leaf: leaf.name, trees.next().get_terminals() )
genomes = re.findall('^(.*?)\.faa$', open( '../../genomes_4_clustering_homologues/tmp/selected.genomes').read(), re.M)

trees = []
for tree in phy.parse('gene_cluster-2.tre', 'newick'):
    trees.append(tree)

branches_df     = pd.DataFrame(index=genomes)
tree_branches   = {}
branch_counter  = 0
tree_counter    = -1
tree_ids        = []
for tree in trees:
    tree_counter += 1
    tree_name     = str(tree_counter)
    tree_branches[tree_name] = []
    tree_ids.append(tree_name)

    tree_leaves = map( lambda leaf: leaf.name , tree.get_terminals() )

    for branch in tree.get_nonterminals():
        if branch.confidence < 0.99:
            continue
        
        branch_counter += 1

        branch_leaves = map( lambda leaf: leaf.name , branch.get_terminals() )

        pb = map( lambda leaf: 1 if leaf in branch_leaves else 0 , tree_leaves )

        branches_df[branch_counter] = pd.Series(pb, tree_leaves)
        tree_branches[tree_name].append(branch_counter)

    if not tree_branches[tree_name]:
        tree_branches.pop(tree_name)
        tree_ids.remove(tree_name)

partitions = [ [] , [] ]
for tree_id in tree_ids:
    partitions[1].append( branches_df[tree_branches[tree_id]].copy() )

def get_incompatible_branches(dfs):
    raw        = []
    normalized = []

    for df1 , df2 in dfs:
    
        intersection = df1.dropna(how='all').index.intersection(df2.dropna(how='all').index)
        incompatible = 0

        df1 = df1.reindex(index=intersection)
        df2 = df2.reindex(index=intersection)
        for branch in df1.columns:
            s1 = set(df1[branch][df1[branch]==1].index)

            for comp_branch in df2.columns:

                s2      = set(df2[comp_branch][df2[comp_branch]==1].index)
                if s1.issuperset(s2) or s1.issubset(s2):
                    continue

                s2_comp = set(df2[comp_branch][df2[comp_branch]==0].index)
                if s1.issuperset(s2_comp) or s1.issubset(s2_comp):
                    continue

                incompatible += 1
                break

        incompatible *= 2.
        raw.append(incompatible)
        normalized.append( incompatible/(2*(len(intersection)-3)) )

    return (raw, normalized)

del(rf_distances)

tree_pairs = list( combinations( partitions[1], 2 ) )

tree_pairs = []
for i in partitions[0]:
    for j in partitions[1]:
        tree_pairs.append([i,j])

num_of_threads = 18
num_of_comparisons = len(tree_pairs)
avg = num_of_comparisons / float(num_of_threads)
datasets = []
last = 0.0
while last < num_of_comparisons:
    datasets.append(tree_pairs[int(last):int(last + avg)])
    last += avg

start_time = time()
pool = Pool(processes=18)
rf_distances = pool.map(get_incompatible_branches, datasets )
pool.close()
pool.join()
print time() - start_time

out = open('rf_distances.cluster2-cluster3', 'wb')
dump(rf_distances, out)
out.close()



plt.figure()
plt.hist(raw3, bins=50, alpha=0.7, histtype='stepfilled', color='y')
plt.hist(raw4, bins=50, alpha=0.7, histtype='stepfilled', color='b')
plt.hist(raw34, bins=50, alpha=0.7, histtype='stepfilled', color='g')
plt.savefig('rf_dists-hist_raw.pdf')
plt.close()

plt.figure()
plt.hist(normal3, bins=50, alpha=0.7, histtype='stepfilled', color='y', normed=True)
plt.hist(normal4, bins=50, alpha=0.3, histtype='stepfilled', color='b', normed=True)
plt.hist(normal34, bins=50, alpha=0.3, histtype='stepfilled', color='g', normed=True)
plt.savefig('rf_dists-hist_normalized_normed-34.pdf', bbox_inches='tight')
plt.close()
