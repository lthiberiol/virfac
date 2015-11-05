from matplotlib import use
use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import re

presence_absence = pd.read_table('../../presence_absence-merged.tab', index_col=0)

clusters = []
for line in open('gene_clusters-0.7-I55').xreadlines():
    clusters.append(line.split())

fig, [ax1, ax2] = plt.subplots( nrows=2 )
ax1.set_title( 'Species per group' )
sns.distplot( presence_absence[clusters[0]].sum(), ax=ax1 )
ax2.set_title( 'Groups per species' )
sns.distplot( presence_absence[clusters[0]].T.sum(), ax=ax2 )

fig.tight_layout()
fig.savefig('cluster0_distribution_test.pdf')

species_per_group = map( 
    lambda position: presence_absence[clusters[position]].sum().values,
    range(50) ) 
tmp = np.asarray( map( 
    lambda position: presence_absence[clusters[position]].T.sum().values / float( len( clusters[position] ) ),
    range(50) ) )
groups_per_species = []
for line in tmp.T:
    groups_per_species.append( line[line > 0] )

fig, [ax1, ax2] = plt.subplots(nrows=2, figsize=(25,15), dpi=300 )

sns.violinplot(species_per_group, ax=ax1, color='Set3')
ax1.set_xticklabels( map( lambda position: 'cluster_%i' %position, range(50) ), rotation='vertical' )
ax1.set_title( 'Species per group' )

sns.violinplot(groups_per_species, ax=ax2, color='Set3')
ax2.set_xticklabels( presence_absence.index, rotation='vertical' )
ax2.set_title( 'Groups per species' )

fig.tight_layout()
fig.savefig('violin_test.pdf')
print 'yeah'


for position in range(50):

    concatenation = {}
    for species in presence_absence.index:
        concatenation[species] = ''

    all_species = set( presence_absence.index )
    for group in clusters[position]:
        print group

        tmp_aln = open('mafft/%s.aln' %group).read()
        tmp_aln = re.sub( '^>\S+\|', '>', tmp_aln, flags=re.M ).split('>')
        tmp_aln.pop(0)

        observed_species = []
        aln_length    = set()
        fucked_up     = False
        for block in tmp_aln:
            block   = block.split('\n') 
            species = block[0]
            seq     = ''.join( block[1:] )

            aln_length.add( len(seq) )

            if species in observed_species:
                fucked_up = True
                break
            else:
                observed_species.append(species)
                concatenation[species] += seq

        if len(aln_length) != 1:
            fucked_up = True
        else:
            aln_length = aln_length.pop()
            for species in all_species.difference( observed_species ):
                concatenation[species] += '-' * aln_length

        if fucked_up:
            print '**WTF!!!'
            break

    out = open('cluster%i_concat.aln' %position, 'wb')
    for species in concatenation:
        if re.match( '-+$', concatenation[species] ):
            continue
        out.write('>%s\n%s\n' %( species, concatenation[species] ) )
    out.close()
