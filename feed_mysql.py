__author__ = 'Thiberio'

import MySQLdb as sql
from os import chdir, listdir
import re

db  = sql.connect(user='root', db='aeromonas')
cursor = db.cursor()

chdir('/Users/Thiberio/work/virulence_factors')

for tmp in listdir('homologous_groups'):
    if not tmp.endswith('.faa'):
        continue

    group_name = tmp.strip('.faa')
    sequences = {'dna':{}, 'aa':{}}

    print group_name

    sql_command = []
    faa = open('homologous_groups/%s' %tmp).read()
    for seq_block in faa.split('>'):
        if seq_block == '':
            continue
        seq_block = seq_block.strip()
        seq_block = seq_block.split('\n')
        header = re.search('^(\S+)\s\[(.+?)\]', seq_block[0], re.M)
        seq = seq_block[1]

        sql_command.append([header.group(1), header.group(2).replace(' ', '_'), seq])
#        sql_command.append('("%s" , "%s" , "%s")' %(header.group(1), header.group(2).replace(' ', '_'), seq))
#    print 'INSERT INTO sequences_test (seq_id , genome , aminoacid) VALUES %s' %' , '.join(sql_command)
    cursor.executemany('INSERT INTO sequences (seq_id , genome , aminoacid) VALUES (%s , %s , %s)',sql_command)
db.commit()



dna_files = listdir('public/DNA')
aa_files = listdir('public/Protein')

for tmp in listdir('public/DNA'):
    if not tmp.endswith('.fna'):
        continue

    genome_name = tmp.strip('.fna')
    if genome_name in sequences:
        print '\t**%s already on db!' %genome_name
        continue
    sequences[genome_name] = {'dna':{}, 'aa':{}}
    print genome_name

    fna = open('public/DNA/%s' %tmp).read()
    for seq_block in fna.split('>'):
        if seq_block == '':
            continue
        seq_block = seq_block.split('\n')
        header = seq_block[0]
        seq_iq = header.strip(genome_name).strip('.')
        seq = ''.join(seq_block[1:]).upper()
        sequences[genome_name]['dna'][seq_iq] = seq
    del(fna)

    if genome_name+'.faa' not in aa_files:
        '**FUDEU, falta correspondência entre arquivos!'
        break
        continue

    faa = open('public/Protein/%s.faa' %genome_name) .read()
    for seq_block in faa.split('>'):
        if seq_block == '':
            continue
        seq_block = seq_block.split('\n')
        header = seq_block[0]
        seq_iq = header.strip(genome_name).strip('.')
        seq = ''.join(seq_block[1:]).upper()
        sequences[genome_name]['aa'][seq_iq] = seq
    del(faa)

    if set(sequences[genome_name]['dna'].keys()) != set(sequences[genome_name]['aa'].keys()):
        print '**Fudeu, o nome das sequências não bate!'
        break

yeah = open('../virulence_factors/genomes_4_clustering_homologues/tmp/selected.genomes').read()
yeah = yeah.split('.faa\n')
pwd