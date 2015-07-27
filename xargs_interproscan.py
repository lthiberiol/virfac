#!/usr/bin/python

from os import system
from sys import argv

out_folder1 = '/work/hgt/bepe/aeromonas/virulence_factors/interproscan'

fasta_folder = '/work/hgt/bepe/aeromonas/virulence_factors/homologous_groups'

fasta = argv[1]+'.faa'

system('interproscan -i ../homologous_groups/%s -appl PfamA-27.0,TIGRFAM-13.0,Panther-9.0,ProDom-2006.1 --iprlookup --goterms -f gff3,html' %fasta)
