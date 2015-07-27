#!/usr/bin/python

from os import system
from sys import argv

out_folder = '/work/hgt/bepe/aeromonas/virulence_factors/profile_annotation/uniref90'

system('hmmsearch -o %s/%s %s /work/databases/blast/nr' %(out_folder, argv[1]+'.hmmout', argv[1]+'.hmm'))
