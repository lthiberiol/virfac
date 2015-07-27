__author__ = 'Thiberio'


from os import system
from sys import argv

out_folder = '~/work/virulence_factors'

system('mafft --auto %s > %s/%s' %(argv[1], out_folder, argv[1].replace('fasta', 'aln')))