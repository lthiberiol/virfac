#!/Users/Thiberio/thiba/bin/python

from os import system
from sys import argv

out_folder = '../trees'

system('fasttree %s > %s/%s' %(argv[1], out_folder, argv[1].replace('.aln', '.tre')))
