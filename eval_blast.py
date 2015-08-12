#!/usr/bin/python
#coding: utf-8
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from statsmodels.sandbox.stats.multicomp import fdrcorrection0, multipletests
from scipy.stats import ttest_ind
import matplotlib.pyplot as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import re
from os import listdir, chdir, system
from Bio.SearchIO import parse
from pickle import load, dump
from sys import exit
import matplotlib.gridspec as gridspec
import scipy.spatial.distance as distance
import scipy.cluster.hierarchy as sch
import multiprocessing
from os.path import isfile

groups = load(open('../homologous_groups-merged.pkl'))

genomes = open('../genomes_4_clustering_homologues/tmp/selected.genomes').read().split('\n')
genomes.pop()
genomes.pop()
for position in range(len(genomes)):
    genomes[position] = genomes[position].strip('.faa')

vfdb_folder = 'vfdb'
control_desc = open('%s/VFDB_T6SS_selected_Protein.headers' %vfdb_folder).read()
control_ids = {}
for element in re.findall('^(\S+)\s(.*?)[,[]', control_desc, re.M):
    control_ids[element[0]] = element[1].strip(' ')

df = pd.read_table('../presence_absence-merged.tab', index_col=0)

#
# evaluate hmmer outputs
#
hmm_positives = {}
for group in groups.keys():

    group = group.replace('&', '-')

    result  = group+'.hmm.hmmout'
    if not isfile('%s/hmm/%s' %(vfdb_folder, result)):
        print '%s not found!' %result
        continue

    result = parse('%s/hmm/%s' %(vfdb_folder, result), 'hmmer3-text').next()

    if not result.hits:
        continue

    best_hit = result.hits[0]
    if best_hit.evalue <= 1e-10:
        hsp = best_hit.hsps[0]
        hmm_positives[group] = {
            'bitscore'  : best_hit.bitscore,
            'bias'      : best_hit.bias,
            'evalue'    : best_hit.evalue,
            'acc'       : hsp.acc_avg,
            'coverage'  : float(hsp.query_end - hsp.query_start) / result.seq_len,
            'desc'      : control_ids[best_hit.id]
        }

out = open('hmm_positives.pkl', 'wb')
dump(hmm_positives, out)
out.close()

vfdb_prot_len = {}
for block in open('vfdb/VFDB_T6SS_selected_Protein.fas').read().split('>'):
    if block == '':
        continue
    block = block.strip()
    vfdb_id = block.split()[0]
    seq = ''.join(block.split('\n')[1:])
    vfdb_prot_len[vfdb_id] = float(len(seq))

def run_blastp(group):
    out_folder1 = 'vfdb/blast'
    out_folder2 = 'control/blast'
    db1 = 'vfdb/VFDB_T6SS_selected_Protein.fas' 
    db2 = 'control/control.faa' 

    fasta_folder = '../homologous_groups'

    fasta = group+'.faa'
    blast_res = group+'.bls'

    system('blastp -query %s/%s -db %s -outfmt 7 -out %s/%s' %(fasta_folder, fasta, db1, out_folder1,blast_res))
    system('blastp -query %s/%s -db %s -outfmt 7 -out %s/%s' %(fasta_folder, fasta, db2, out_folder2,blast_res))

pool = multiprocessing.Pool(processes=15)
pool.map(run_blastp, hmm_positives.keys())
pool.close()
pool.join()

all = {}
false_positives = []
true_positives = {}
test_pvalue = {}
control_folder = 'control/blast'
for hg_id in hmm_positives:
    
    vfdb_blast = open('%s/blast/%s.bls' %(vfdb_folder, hg_id)).readlines()
    vfdb_parameters =  {'desc':[], 'scores':[], 'identity':[], 'coverage':[], 'hit coverage':[]}
    for line in vfdb_blast:
        if line.startswith('#'):
            flag = True
            regex = re.search('\[(.*)\]\s\|.*\((\d+)\)$', line, re.M)
            if regex:
                query_length = float(regex.group(2))
            continue
        if flag:
            flag        = False
            line        = line.split()
            hit_id      = line[1]
            bitscore    = float(line[-1])
            identity    = float(line[2])
            aln_length  = int(line[7]) - int(line[6])
            coverage    = aln_length/query_length
            hit_aln_length  = int(line[9]) - int(line[8])
            hit_coverage    = hit_aln_length/vfdb_prot_len[hit_id]
            vfdb_parameters['desc'        ].append(control_ids[hit_id])
            vfdb_parameters['scores'      ].append(bitscore)
            vfdb_parameters['identity'    ].append(identity)
            vfdb_parameters['coverage'    ].append(coverage)
            vfdb_parameters['hit coverage'].append(hit_coverage)

    control_blast = open('%s/%s.bls' %(control_folder, hg_id)).readlines()
    control_parameters = {'hits' : [], 'scores' : [], 'identity':[], 'coverage':[]}
    for line in control_blast:
        if line.startswith('#'):
            flag = True
            regex = re.search('\[(.*)\]\s\|.*\((\d+)\)$', line, re.M)
            if regex:
                query_length = float(regex.group(2))
            continue
        if flag:
            flag        = False
            line        = line.split()
            hit_id      = line[1]
            bitscore    = float(line[-1])
            identity    = float(line[2])
            aln_length  = int(line[7]) - int(line[6])
            coverage    = aln_length/query_length
            control_parameters['hits'    ].append(hit_id)
            control_parameters['scores'  ].append(bitscore)
            control_parameters['identity'].append(identity)
            control_parameters['coverage'].append(coverage)

    all[hg_id] = {'vfdb':vfdb_parameters, 'control':control_parameters}

    if not control_parameters['hits']:
        true_positives[hg_id] = vfdb_parameters
        continue
    
    if len(vfdb_parameters['desc']) < 3 or len(control_parameters['hits']) < 3:
        if np.min(vfdb_parameters['scores']) > 1.5*np.max(control_parameters['scores']):
            true_positives[hg_id] = vfdb_parameters
        else:
            false_positives.append(hg_id)
    else:
        (t, p) = ttest_ind(vfdb_parameters['scores'], control_parameters['scores'])
        pvalue = p/2 #one tailed
        test_pvalue[hg_id]          = vfdb_parameters
        test_pvalue[hg_id]['pvalue']= pvalue
        if t < 0:
            false_positives.append(hg_id)

tmp_pvalue      = []
sorted_hg_ids   = []
for hg_id in test_pvalue:
    tmp_pvalue.append(test_pvalue[hg_id]['pvalue'])
    sorted_hg_ids.append(hg_id)
    test_pvalue[hg_id].pop('pvalue')

corrected_pvalue= fdrcorrection0(tmp_pvalue)[1]
for position in range(len(corrected_pvalue)):
    hg_id = sorted_hg_ids[position]
    if corrected_pvalue[position] <= 0.05 and hg_id not in false_positives:
        true_positives[hg_id] = test_pvalue[hg_id]
    elif corrected_pvalue[position] > 0.05 and hg_id not in false_positives:
        false_positives.append(hg_id)

out = open('true_positives.pkl', 'wb')
dump(true_positives, out)
out.close()

agreement = []
groups_description = {}
descriptions = set()
for hg_id in true_positives:
    blank_spaces = ' '*(30-len(hg_id))
    count = float(true_positives[hg_id]['desc'].count(hmm_positives[hg_id]['desc']))
    agreement.append(count/len(true_positives[hg_id]['desc']))
    descriptions.add(hmm_positives[hg_id]['desc'])
    print '%s:%s %s (%.2f%s)' %(hg_id, blank_spaces, hmm_positives[hg_id]['desc'], 100 * count/len(true_positives[hg_id]['desc']), '%')

df_hg = pd.DataFrame(index=genomes, columns=true_positives.keys())
for hg_id in true_positives:
    for genome in groups[hg_id]:
        df_hg.loc[genome][hg_id] = []
        for protein in groups[hg_id][genome]:
            df_hg.loc[genome][hg_id].append(protein)

df_binary = pd.DataFrame(index=genomes, columns=true_positives.keys(), data=0)
for hg_id in true_positives:
    for genome in groups[hg_id]:
        df_binary.loc[genome][hg_id] = 1

all_identities      = []
all_coverage        = [] 
positive_identities = []
positive_coverages  = []
negative_identities = []
negative_coverages  = []
scores = {'vfdb':[], 'control':[]}
positive_scores = {'vfdb':[], 'control':[]}
negative_scores = {'vfdb':[], 'control':[]}
for n in all:
    identity = all[n]['vfdb']['identity']
    coverage = all[n]['vfdb']['coverage']
    if np.mean(coverage) > 1:
        continue
    tmp1 = np.mean(all[n]['vfdb']['scores'])
    tmp2 = 0.0
    if all[n]['control']['scores']:
        tmp2 = np.mean(all[n]['control']['scores'])
    scores['vfdb'   ].append(tmp1)
    scores['control'].append(tmp2)
    if n in true_positives:
        positive_scores['vfdb'].append(tmp1)
        positive_scores['control'].append(tmp2)
        positive_identities.append(np.mean(identity))
        positive_coverages.append(np.mean(coverage))
    else:
        negative_scores['vfdb'].append(tmp1)
        negative_scores['control'].append(tmp2)
        negative_identities.append(np.mean(identity))
        negative_coverages.append(np.mean(coverage))
    all_identities.append(np.mean(identity))
    all_coverage.append(  np.mean(coverage))


################################################################################
########## print bistscore and coverage histogram of each hg best hits on VFDB #
################################################################################
# define figure size and dpi
plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] = (15, 15)
# set X and Y axis
fig, axScatter = pl.subplots()
fig.suptitle("Average BLAST bitscores", fontsize=30)
axScatter.set_xlabel('T6SS')
axScatter.set_ylabel('Control')
axScatter.scatter(positive_scores['vfdb'], positive_scores['control'], alpha=0.5, c='b', edgecolor='none', s=50, label='T3SS better scores')
axScatter.scatter(negative_scores['vfdb'], negative_scores['control'], alpha=0.5, c='k', edgecolor='none', s=50, label='Control better scores')
axScatter.grid()
axScatter.legend(scatterpoints=1)
divider = make_axes_locatable(axScatter)
axHistx=divider.append_axes('top', 1.2, pad=0.1, sharex=axScatter)
axHisty=divider.append_axes('right', 1.2, pad=0.1, sharey=axScatter)
pl.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(), visible=False)
axHistx.hist(scores['vfdb'], bins=100, color='k', alpha=0.7, edgecolor='white')
axHisty.hist(scores['control'], bins=100, orientation='horizontal', color='k', alpha=0.7, edgecolor='white')
axHisty.set_ylim(bottom=0)
axHistx.set_yticks([])
axHisty.set_xticks([])
fig.tight_layout()
fig.savefig('blast_scores.pdf')
################################################################################

################################################################################
########## print bistscore and coverage histogram of each hg best hits on VFDB #
################################################################################
# set X and Y axis
x = all_identities
y = all_coverage

fig, axScatter = plt.subplots()
fig.suptitle("Mean identities and coverage of putative T3SSs", fontsize=30)
axScatter.set_xlabel('Identity (%)')
axScatter.set_ylabel('Coverage (% of query)')
#axScatter.scatter(x, y, alpha=0.25, c='k', edgecolor='none', s=50)
axScatter.scatter(positive_identities, positive_coverages, alpha=0.5, c='b', edgecolor='none', s=50)
#axScatter.scatter(negative_identities, negative_coverages, alpha=0.5, c='k', edgecolor='none', s=50)
axScatter.grid()
divider = make_axes_locatable(axScatter)
axHistx=divider.append_axes('top', 1.2, pad=0.1, sharex=axScatter)
axHisty=divider.append_axes('right', 1.2, pad=0.1, sharey=axScatter)
plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(), visible=False)
axHistx.hist(x, bins=100, color='k', alpha=0.7, edgecolor='white')
axHisty.hist(y, bins=100, orientation='horizontal', color='k', alpha=0.7, edgecolor='white')
axHisty.set_ylim(bottom=0)
axHistx.set_xlim(right=100)
axHistx.set_yticks([])
axHisty.set_xticks([])
fig.tight_layout()
fig.savefig('mean_identity_coverage-no_black.pdf')
################################################################################
