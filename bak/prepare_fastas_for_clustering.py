#!/usr/bin/python
#coding: utf-8

from os import listdir
import re

source = '/work/hgt/bepe/aeromonas/original_files_gave_me_by_matt/Public_Genomes/'
output_folder = '/work/hgt/bepe/aeromonas/virulence_factors/new_genomes/'

all_faas = open('new_genomes-verified_to_use.list').read().split('\n')

for faa in all_faas:
    if faa.startswith('#'):
        continue

    print faa
    
    entry  = open(source+faa+'.faa').readlines()
    output = open(output_folder+faa+'.faa', 'wb')
    taxon_name = faa.replace('_', ' ')

    for line in entry:
        if not line.startswith('>'):
            output.write(line)
            continue

        line = line.strip('\n')
        protein_id = re.search('[.|](.*)$', line, re.M)

        output.write('>%s [%s]\n' %(protein_id.group(1), taxon_name))
    
    output.close()
