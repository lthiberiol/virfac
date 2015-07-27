#!/Users/Thiberio/thiba/bin/python
from sys import argv, exit
from  os import system

if not argv[1].endswith('.phy'):
    exit()

system('puzzle %s < ../puzzle_params' %argv[1])
