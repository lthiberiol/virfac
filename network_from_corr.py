__author__ = 'Thiberio'

from matplotlib import use
use('Agg')
from matplotlib import pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
from os import chdir, listdir, system
from itertools import combinations
from scipy.spatial.distance import squareform
from pickle import load
from math import sqrt
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.figsize'] = (15, 15)

chdir('/Volumes/Macintosh HD 2/thiberio/virulence_factors')

groups = load( open( 'homologous_groups-merged.pkl' ) )
corr_df = pd.read_table('significant_group_correlations-based_on_distances.tab', index_col=0)
indexes = corr_df.index[:1000]
smaller_df = corr_df.loc[indexes][indexes].copy()

writable_df = corr_df.copy()
writable_df[writable_df < 0] = 0
writable_df.to_csv('mcl_input.abc', sep='\t')
system('mcxarray -data mcl_input.abc -co 0.9 -skipr 1 -skipc 1 -o mcl_input-0.9.mci -write-tab mcl_input.dict --spearman -t 10')
system('mcx query -imx mcl_input.mci -o mcxq.out --vary-correlation -t 10')

graph = nx.Graph()

graph.add_nodes_from(corr_df.index)

weights = squareform(corr_df.fillna(0).values)
edges_between = list( combinations( corr_df.index, 2))
for n in range(len(weights)):
    edge_weight = weights[n]
    if edge_weight >= 0.7:
        node1 = edges_between[n][0]
        node2 = edges_between[n][1]

        graph.add_edge(node1, node2, weight = edge_weight)

nx.graph_number_of_cliques(graph)
print nx.info(graph)
out = open('graph.dimacs', 'wb')
for line in graph.degree(weight='weight').iteritems():
    out.write('n %s %f\n' %line)
for line in graph.edges():
    out.write('e %s %s\n' %line)

cliques = nx.find_cliques(graph)
cliques = list(cliques)

a = smaller_df.as_matrix()

color_intensity = []
node_sizes      = []
graph_degree = graph.degree(weight='weight')
for node in graph.nodes():
    color_intensity.append( graph_degree[node] )
    node_sizes.append( len( groups[node] ) )

print '\t**Generating spring layout ...'
graph_layout = nx.spring_layout(graph, iterations=50)
graph_layout2 = nx.spectral_layout(graph)
print '... done!'

plt.figure()
plt.xticks([])
plt.yticks([])

print '\t**Drawing nodes ...'
nx.draw_networkx_nodes(graph, pos=graph_layout2, node_color=color_intensity, with_labels=False, alpha=0.75, node_size= 100, cmap=plt.get_cmap('Blues'))
print '\t**Saving pdf (no edges) ...'
plt.savefig('correlation_network5-no_edges.pdf', bbox_inches='tight')
print '\t**Drawing edges ...'
nx.draw_networkx_edges(graph, pos=graph_layout, with_labels=False, alpha=0.3)
print '\t**Saving final pdf ...'
plt.savefig('correlation_network3.pdf', bbox_inches='tight')
print '... done!'
plt.close()

a=nx.find_cliques(graph)
nx.graph_clique_number(graph)
nx.graph_number_of_cliques(graph)
nx.density(graph)
print nx.info(graph)
nx.average_clustering(graph, weight='yeah')
nx.get_edge_attributes(graph, 'weights')