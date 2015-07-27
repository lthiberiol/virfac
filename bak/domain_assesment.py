from matplotlib import use
use('Agg')
import pandas as pd
import pylab as pl
import numpy as np
import re
from pickle import load, dump

entry = open('homologous_groups.list').read()
groups = load(open('homologous_groups.pkl'))
t3ss = load(open('true_positives.pkl'))

go_groups = {}
proteins = []
go_terms = set()
for group in t3ss:
    gff = open('interproscan/%s.faa.gff3' %group).read()
    gff_frags = re.findall('(#sequence-region\s(\S+)(?:.|\n)*?)\n#', gff)
    if gff_frags:
        go_groups[group] = {}
    for fragment in gff_frags:
        prot_name = fragment[1]
        go = set(re.findall('"(GO:\d+)"', fragment[0]))
        if go:
            go_terms.update(go)
            proteins.append(prot_name)
            go_groups[group][prot_name] = go
    if group in go_groups and not go_groups[group]:
        go_groups.pop(group)

go_count = {}
for go in go_terms:
    go_count[go] = 0.0

go_terms_present ={}
for group in go_groups:
    go_terms_present[group] = set()
    for protein_go in go_groups[group].values():
        go_terms_present[group].update(protein_go)
    for go in go_terms:
        if go in  go_terms_present[group]:
            go_count[go] += 1

entry = open('go_terms.desc').read()
entry = entry.split('\n')
go_desc = {}
for line in entry:
    desc, term_type, acc = line.split('\t')
    if term_type not in go_desc:
        go_desc[term_type] = {}
    go_desc[term_type][acc] = desc

for term_type in go_desc:
    pl.figure(figsize=(20,20), dpi=300)
    ticks = []
    values = []
    for  go in go_desc[term_type]:
        values.append(go_count[go])
        ticks.append(go_desc[term_type][go])
    pl.bar(np.arange(len(ticks)), values)
    pl.xticks(np.arange(len(ticks)), ticks, rotation=45, fontsize=15, ha='right')
    pl.tight_layout()
    pl.savefig('ipad_sync/%s-bar.png' %term_type)

