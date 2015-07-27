#!/Users/Thiberio/thiba/bin/python
#coding: utf-8

from pickle import load
from sys import argv
import re

puzzle = open(argv[1]).read()

ids = re.findall('^(_\d+_)', puzzle, re.M)
codes = load(open('../mafft-seq_codes/%s' %argv[1].replace('.phy.dist', '-seq_codes.pkl')))
inv_codes = {v: k for k, v in codes.items()}

species = []
for identifier in ids:
    species.append(inv_codes[identifier])

puzzle = puzzle.split('\n')
puzzle.pop(0)
puzzle.pop()

distances = {}
for line in puzzle:
    search_actual = re.search('^(_\d+_)', line, re.M)

    if search_actual:
        actual = search_actual.group(1)
        distances[inv_codes[actual]] = []
        line = line.split()[1:]
    else:
        line = line.split()
    
    for element in line:
        distances[inv_codes[actual]].append(element)

out = open('../distances/%s.tab' %argv[1], 'wb')
out.write('gene_names\t%s\n' %'\t'.join(species))
for species_name in species:
    out.write('%s\t%s\n' %(species_name, '\t'.join(distances[species_name])))
out.close()
