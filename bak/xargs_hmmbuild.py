#!/usr/bin/python

from os import system
from sys import argv

out_folder = '/work/hgt/bepe/aeromonas/virulence_factors/hmm_profiles'

system('hmmbuild --amino %s/%s %s' %(out_folder, argv[1].replace('.aln', '.hmm'), argv[1]))
