__author__ = 'Thiberio'


from pickle import load, dump
import MySQLdb as sql
import pandas as pd
from Bio.SeqIO import convert
from os import listdir, chdir
import re
from os.path import getsize
from scipy.stats import spearmanr, pearsonr
from scipy.spatial.distance import squareform
from itertools import combinations
from progressbar import ProgressBar
from copy import deepcopy
pbar = ProgressBar()

chdir('/Volumes/Macintosh HD 2/thiberio/virulence_factors/')

db  = sql.connect(user='root', db='aeromonas')
cursor = db.cursor()

groups = load(open('homologous_groups-merged.pkl'))

###############################################################################
# create #
###############################################################################
pb = pd.read_table('presence_absence-merged.tab', index_col=0, nrows=105)

for group in df.sum()[df.sum() >= 10].index:
    sql_query = []
    for genome in groups[group]:
        for protein in groups[group][genome]:
            sql_query.append('(genome="%s" and seq_id="%s")' %(genome, protein))
    cursor.execute('select seq_id, genome, aminoacid from sequences where %s' %' or '.join(sql_query))

    out = open('fastas/%s.faa' %group.replace('&', '-'), 'wb')
    for query_result in cursor.fetchall():
        out.write('>%s|%s\n%s\n' %query_result)
    out.close()

for mafft_output in listdir('mafft/'):
    entry = open('mafft/%s' %mafft_output).read()
    headers = re.findall('^>(.*)$', entry, re.M)
    counter = 0
    codes = {}
    for header in headers:
        counter += 1
        entry = entry.replace(header, '_%s_' %str(counter))
        codes[header] = '_%s_' %str(counter)

    out = open('mafft_abbr/%s' %mafft_output, 'wb')
    out.write(entry)
    out.close()

    out = open('mafft-seq_codes/%s' %mafft_output.replace('.aln', '-seq_codes.pkl'), 'wb')
    dump(codes, out)
    out.close()

phylip_files = listdir('phylip')
for mafft_output in listdir('mafft_abbr/'):
    if mafft_output.replace('.aln', '.phy') in phylip_files:
        continue
    if getsize('mafft_abbr/%s' %mafft_output):
        convert('mafft_abbr/%s' %mafft_output, 'fasta', 'phylip/%s' %mafft_output.replace('.aln', '.phy'), 'phylip')


######################################################
########################### finally, assess distances!
######################################################
count = 0
max_count = len(listdir('tree_puzzle/distances'))
nice_families = {}
usefull_genes = {}
pbar = ProgressBar()
for dist_table in pbar(listdir('tree_puzzle/distances')):

    flag = False

    if not dist_table.endswith('.tab'):
        continue

    group = dist_table.replace('.phy.dist.tab', '')
    group = group.replace('-', '&')

    df = pd.read_table('tree_puzzle/distances/%s' %dist_table, index_col=0)

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

    # usefull_genes[group] = []
    # for genome in groups[group]:
    #     if len(groups[group][genome]) == 1:
    #         usefull_genes[group].append( '%s|%s' %(groups[group][genome][0], genome) )
    #         continue
    #
    #     minimum_dist_paralog = {'gene_name': '', 'dist_sum': 1000}
    #
    #     for protein in groups[group][genome]:
    #         tmp_dists = df.loc['%s|%s' %(protein, genome)]
    #
    #         tmp_dist_sum = tmp_dists.sum()
    #         if tmp_dist_sum < minimum_dist_paralog['dist_sum']:
    #             minimum_dist_paralog['dist_sum']  = tmp_dist_sum
    #             minimum_dist_paralog['gene_name'] = '%s|%s' %(protein, genome)
    #
    #         non_null_dists = tmp_dists[tmp_dists > 0]
    #         non_null_dists.sort()
    #         for index in non_null_dists.index:
    #             if index.endswith('genome'):
    #                 continue
    #             break
    #
    #         for paralog in groups[group][genome]:
    #             if protein == paralog:
    #                 continue
    #
    #             if tmp_dists['%s|%s' %(paralog, genome)] > tmp_dists[index]:
    #                 #print 'Paralogs do not form a tight cluster!'
    #                 flag = True
    #                 break
    #
    #         if flag:
    #             break
    #
    #         usefull_genes[group].append( minimum_dist_paralog['gene_name'] )
    #
    #     if flag:
    #         usefull_genes.pop(group)
    #         break
    #
    # if not flag:
    #     nice_families.append(group)

species_combinations = []
for pair in combinations(pb.index, 2):
    species_combinations.append( frozenset(pair))
all_species_distances = pd.DataFrame(index=nice_families.keys(), columns=species_combinations)

pbar = ProgressBar()
for group in pbar(nice_families.keys()):

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
all_species_distances.to_csv('all_all_species_distances.tab', sep='\t')
out = open('all_species_distances.pkl', 'wb')
dump(all_species_distances, out)
out.close()
'--'.join(all_species_distances.columns[0])

group_combinations  = list( combinations(nice_families.keys(), 2) )
condensed_corr = []
condensed_pval = []
pbar = ProgressBar()
for pair in pbar(group_combinations):
    (rho, p) = spearmanr( all_species_distances.loc[list( pair )].dropna(axis=1, how='any').T )
    condensed_corr.append(rho)
    condensed_pval.append(p)


#distances_corr = all_species_distances.T.corr(method='spearman')
#distances_corr.to_csv('gene_family_distances_correlations.tab', sep='\t')
distances_corr = pd.read_table('gene_family_distances_correlations.tab', index_col=0)

from statsmodels.sandbox.stats.multicomp import MultiComparison
a = all_species_distances.loc[all_species_distances.index[:5]][all_species_distances.columns[:5]]
b=a.T.corr(method='pearson')