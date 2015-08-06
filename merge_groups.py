from os import chdir, listdir, system
from os.path import isfile, getsize
from pickle import load,dump
from Bio import Phylo as phy
import pandas as pd
from copy import deepcopy
import multiprocessing
from Bio.SearchIO import parse
from itertools import combinations

chdir('/Volumes/Macintosh HD 2/thiberio/virulence_factors')
hmm_positives = load(open('hmm_positives.pkl'))
true_positives = load(open('true_positives.pkl'))
groups = load(open('../homologous_groups.pkl'))
df = pd.read_table('presence_absence.tab', index_col=0)
hm21 = 'A_veronii_Hm21'

same_desc = {}
for group in true_positives:
    desc = hmm_positives[group]['desc']
    desc = desc.replace('/', '-').replace(' ', '_')
    if desc not in same_desc:
        same_desc[desc] = []
    same_desc[desc].append(group)

for desc in same_desc:

    if len(same_desc[desc]) == 1:
        continue

    out = open('group_merge/%s.faa' %desc.replace('/', '-').replace(' ', '_'), 'wb')
    for group in same_desc[desc]:
        entry = open('../homologous_groups/%s.faa' %group).read()
        out.write(entry.replace('>', '>%s|' %group).strip('\n'))
        out.write('\n')
    out.close()

################################################################################
############### run HMMSEARCH for all groups that should be tested for merging #
################################################################################
def run_hmm(group):
    system('hmmsearch -o profile_annotation/%s.hmmout ../hmm_profiles/%s.hmm ../../blast_dbs/uniref90.fasta' %(group, group))

pool = multiprocessing.Pool(processes=15)
pool.map(run_hmm, true_positives.keys())
pool.close()
pool.join()
################################################################################
uniref_homologues = {}
for group in true_positives.keys():
    print group

    hmmsearch = parse('profile_annotation/%s.hmmout' %group, 'hmmer3-text').next()

    if not hmmsearch.hits:
        continue
    
    best_hit = hmmsearch.hits[0]
    bitscore_threshold = best_hit.bitscore * 0.7

    uniref_homologues[group] = set()
    for hit in hmmsearch.hits:
        if hit.evalue > 1e-15 or hit.bitscore < bitscore_threshold:
            break
        else:
            uniref_homologues[group].add(hit.id)

    if not uniref_homologues[group]:
        uniref_homologues.pop(group)

yeah = []
for pair in combinations(uniref_homologues.keys(), 2):
    if uniref_homologues[pair[0]].intersection(uniref_homologues[pair[1]]):
        flag = True
        pair = set(pair)
        for pos in range(len(yeah)):
            if yeah[pos].intersection(pair):
                yeah[pos].update(pair)
                flag = False
                break
        if flag:
            yeah.append(pair)

def chunks(entry, delim='\n'):
    buf = bytearray(), 
    while True:
        c = entry.read(1)
        if c == '': 
            yield str(buf)
            return
        print c
        buf.append(c)
        if c == delim: 
            yield str(buf)
            buf = bytearray()

print 'yeah'
best_hit.bitscore

################################################################################
###### run MSA and tree for all descriptions that should be tested for merging #
################################################################################
def run_tree(fasta):
    if isfile('group_merge/%s.faa' %fasta):
        system('mafft --auto --reorder group_merge/%s.faa > group_merge/mafft/%s.aln' %(fasta, fasta))
        system('fasttree -wag -gamma group_merge/mafft/%s.aln > group_merge/trees/%s.tre' %(fasta, fasta))

pool = multiprocessing.Pool(processes=15)
pool.map(run_tree, same_desc.keys())
pool.close()
pool.join()
################################################################################

should_merge = {}
# should_merge = []
for tree_name in listdir('group_merge/trees'):

    if not getsize('group_merge/trees/%s' %tree_name):
        continue

    desc = tree_name.split('.tre')[0]
    tree = phy.parse('group_merge/trees/%s' %tree_name, 'newick')
    tree = tree.next()

#    if tree.total_branch_length() > 0:
#        tree.root_at_midpoint()

    tmp_groups = {}
    for leaf in tree.get_terminals():
        group_name = leaf.name.split('|')[0]
        if group_name not in tmp_groups:
            tmp_groups[group_name] = []
        tmp_groups[group_name].append(leaf)

    if len(tmp_groups.keys()) == 1:
        continue

    tree.root_at_midpoint()

    for group in tmp_groups:
        if not tree.is_monophyletic(tmp_groups[group]):
#            should_merge.append(desc)
#            break
            common_ancestor = tree.common_ancestor(tmp_groups[group])
            tmp_should_merge = set()

            for leaf in common_ancestor.get_terminals():
                group_name = leaf.name.split('|')[0]
                tmp_should_merge.add(group_name)
            # tmp_should_merge = frozenset(tmp_should_merge)

            if desc not in should_merge:
                should_merge[desc] = set()

            to_remove = set()
            for mergeable in should_merge[desc]:
                if mergeable in to_remove:
                    continue

                if tmp_should_merge.intersection(mergeable):
                    tmp_should_merge.update(mergeable)
                    if mergeable != tmp_should_merge:
                        to_remove.add(mergeable)

            should_merge[desc].add(frozenset(tmp_should_merge))
            for removable in to_remove:
                should_merge[desc].remove(removable)

all_mergeable_groups = set()
for desc in should_merge:
    for mergeable in should_merge[desc]:
        all_mergeable_groups.update(mergeable)

new_groups  = {}
t3ss_groups = []
for group in groups:
    if group not in all_mergeable_groups:
        new_groups[group] = deepcopy(groups[group])
        if group in true_positives:
            t3ss_groups.append(group)

for desc in should_merge:
    for mergeable in should_merge[desc]:
        new_group_name = '&'.join(mergeable)
        new_groups[new_group_name] = {}
        t3ss_groups.append(new_group_name)
        for group in mergeable:
            if group in new_groups:
                print 'porréssa?!?!'
                break
            for species in groups[group]:
                if species not in new_groups[new_group_name]:
                    new_groups[new_group_name][species] = deepcopy(groups[group][species])
                else:
                    new_groups[new_group_name][species].extend(groups[group][species])

out = open('homologous_groups-merged.pkl', 'wb')
dump(new_groups, out)
out.close()
out = open('t3ss_groups-merged.pkl', 'wb')
dump(t3ss_groups, out)
out.close()

new_df = pd.DataFrame(index=df.index, columns=new_groups.keys(), data=0)
for column in new_df.columns:
    for index in new_df.index:
        if index in new_groups[column]:
            new_df.loc[index][column] = 1
out = open('presence_absence-merged.tab', 'wb')
new_df.to_csv('presence_absence-merged.tab', sep='\t')

mergeable_units_ref = []
mergeable_units     = []
out = open('group_merge/paralogs-mid_point_rooted', 'wb')
for desc in should_merge:
    out.write('%s (%i)\n' %(desc, len(same_desc[desc])))

    for mergeable in should_merge[desc]:
        mergeable = list(mergeable)
        mergeable_units.append(mergeable)

        hm21_flag = False
        if df.loc[hm21][mergeable].sum() == 1:
            mergeable_units_ref.append(mergeable)
            hm21_flag = True

        out.write('\t%s' %', '.join(mergeable))
        if hm21_flag:
            out.write(' (Hm21 reference)\n')
        else:
            out.write('\n')

        df_sum = df[mergeable].T.sum()
        for genome in df_sum[df_sum > 1].index:
            out.write('\t\t%s (%i)\n' %(genome, df_sum[genome]))
        out.write('\n')
out.close()

#
# aexT test!
#
# reference (Hm21)
aext_groups = df[same_desc['aexT']]
reference_aext = aext_groups.loc['A_veronii_Hm21'][aext_groups.loc['A_veronii_Hm21'] == 1].index[0]
hm21_aext = groups[reference_aext]['A_veronii_Hm21'][0]

cursor.execute('SELECT seq_id, genome, aminoacid FROM sequences where seq_id="%s" and genome="%s";' %(hm21_aext, 'A_veronii_Hm21'))
query_result = cursor.fetchone()

out = open('group_merge/referece_aext.faa', 'wb')
out.write('>%s|%s\n%s' %query_result)
out.close()

chdir('group_merge')
genomes_with_paralogs = df.loc[df[same_desc['aexT']].T.sum() > 1].index
for genome in genomes_with_paralogs:
    out = open('%s_aexT.faa' %genome, 'wb')
    for group in same_desc['aexT']:
        if genome not in groups[group]:
            continue
        for protein in groups[group][genome]:
            cursor.execute('SELECT seq_id, aminoacid FROM sequences where seq_id="%s" and genome="%s";' %(protein, genome))
            out.write('>%s\n%s\n' %cursor.fetchone())
    out.close()

    system('mummer -maxmatch reference_aexT.faa %s_aexT.faa > %s.mums' %(genome, genome))
    system('mummerplot --medium --png -c --prefix=%s %s.mums' %(genome, genome))


#
# hey ho! everybody now!
#
for desc in should_merge_ref:
    desc_groups = df[same_desc[desc]]
    reference_desc = desc_groups.loc[hm21][desc_groups.loc[hm21] == 1].index[0]
    hm21_desc = groups[reference_desc][hm21][0]
    cursor.execute('SELECT seq_id, aminoacid FROM sequences where seq_id="%s" and genome="%s";' %(hm21_desc, hm21))

    out = open('reference_%s.faa' %desc, 'wb')
    out.write('>%s\n%s' %cursor.fetchone())
    out.close()

    genomes_with_paralogs = df.loc[df[same_desc[desc]].T.sum() > 1].index
    for genome in genomes_with_paralogs:
        out = open('%s_%s.faa' %(genome, desc), 'wb')
        for group in same_desc[desc]:
            if genome not in groups[group]:
                continue
            for protein in groups[group][genome]:
                cursor.execute('SELECT seq_id, aminoacid FROM sequences where seq_id="%s" and genome="%s";' %(protein, genome))
                out.write('>%s\n%s\n' %cursor.fetchone())
        out.close()

        system('mummer -maxmatch reference_%s.faa %s_%s.faa > %s_%s.mums' %(desc, genome, desc, genome, desc))
        system('mummerplot --medium --png -c --prefix=%s_%s %s_%s.mums' %(genome, desc, genome, desc))

yeah = []
#chdir('group_merge/mergeable_groups_hm21_aln')
for mergeable in mergeable_units_ref:
    desc_groups = df[mergeable]
    reference_desc = desc_groups.loc[hm21][desc_groups.loc[hm21] == 1].index[0]
    hm21_desc = groups[reference_desc][hm21][0]
    cursor.execute('SELECT seq_id, aminoacid FROM sequences where seq_id="%s" and genome="%s";' %(hm21_desc, hm21))

    out = open('reference_%s.faa' %'-'.join(mergeable), 'wb')
    out.write('>%s\n%s' %cursor.fetchone())
    out.close()

    genomes_with_paralogs = df.loc[df[mergeable].T.sum() > 1].index
    for genome in genomes_with_paralogs:
        yeah.append('%s_%s.faa' %(genome, '-'.join(mergeable)))
        out = open('%s_%s.faa' %(genome, '-'.join(mergeable)), 'wb')
        for group in mergeable:
            if genome not in groups[group]:
                continue
            for protein in groups[group][genome]:
                cursor.execute('SELECT seq_id, aminoacid FROM sequences where seq_id="%s" and genome="%s";' %(protein, genome))
                out.write('>%s\n%s\n' %cursor.fetchone())
        out.close()

        system('mummer -maxmatch reference_%s.faa %s_%s.faa > %s_%s.mums' %('-'.join(mergeable), genome, '-'.join(mergeable), genome, '-'.join(mergeable)))
        system('mummerplot --medium --png --prefix=%s_%s %s_%s.mums' %(genome, '-'.join(mergeable), genome, '-'.join(mergeable)))
