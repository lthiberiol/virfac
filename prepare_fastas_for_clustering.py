#!/usr/bin/python
#coding: utf-8

from os import listdir
import re

source = '/Volumes/Macintosh HD 2/thiberio/aeromonas_genomes/DNA_ORFs/'
output_folder = '/Volumes/Macintosh HD 2/thiberio/aeromonas_genomes/DNA_ORFs/formated/'

all_faas = listdir('genomes_4_clustering')

count = 0
for faa in all_faas:
    faa = faa.replace('.faa', '')

    count += 1
    print count, faa
    
    entry  = open(source+faa+'.fna').readlines()
    output = open(output_folder+faa+'.fna', 'wb')
    taxon_name = faa.replace('_', ' ')

    for line in entry:
        if not line.startswith('>'):
            output.write(line)
            continue

        line = line.strip('\n')
        protein_id = re.search('[.|](.*)$', line, re.M)

        output.write('>%s [%s]\n' %(protein_id.group(1), taxon_name))
    
    output.close()
