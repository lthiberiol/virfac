#!/usr/bin/python

from os import system
from sys import argv

out_folder = '../alignments'

system('mafft --auto %s > %s/%s' %(argv[1], out_folder, argv[1].replace('.faa', '.aln')))
