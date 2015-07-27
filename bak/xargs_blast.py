#!/usr/bin/python

from os import system
from sys import argv

out_folder1 = '/work/hgt/bepe/aeromonas/virulence_factors/vfdb/blast'
out_folder2 = '/work/hgt/bepe/aeromonas/virulence_factors/control/blast'
db1 = '/work/hgt/bepe/aeromonas/virulence_factors/vfdb/VFDB_t3ss_Protein.fas' 
db2 = '/work/hgt/bepe/aeromonas/virulence_factors/control/control.faa' 

fasta_folder = '/work/hgt/bepe/aeromonas/virulence_factors/homologous_groups'
fasta = argv[1].replace('.hmmout', '.faa')
blast_res = argv[1].replace('.hmmout', '.bls')

system('blastp -query %s/%s -db %s -outfmt 7 -out %s/%s' %(fasta_folder, fasta, db1, out_folder1,blast_res))
system('blastp -query %s/%s -db %s -outfmt 7 -out %s/%s' %(fasta_folder, fasta, db2, out_folder2,blast_res))
