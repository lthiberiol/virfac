from matplotlib import  use
use('Agg')
import re
import pylab as pl
import numpy as np

entry = open('homologous_groups.list').read()

size = re.findall('size=(\d+)\staxa=\d+\s', entry)
taxa = re.findall('size=\d+\staxa=(\d+)\s', entry)

for n in range(len(size)):
    size[n] = int(size[n])
    taxa[n] = int(taxa[n])

pl.figure(figsize=(15,15),dpi=300)
pl.hist(size, 100)
pl.tight_layout()
pl.savefig('ipad_sync/size_mcl.pdf')

pl.figure(figsize=(15,15),dpi=300)
pl.hist(taxa, 105)
pl.tight_layout()
pl.savefig('ipad_sync/taxa_mcl.pdf')
