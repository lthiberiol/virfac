__author__ = 'Thiberio'

import MySQLdb as sql
from os import chdir, listdir
import re

db  = sql.connect(user='root', db='aeromonas')
cursor = db.cursor()

cursor.execute('''
CREATE TABLE `sequences_all` (
`id` int(11) NOT NULL AUTO_INCREMENT,
`seq_id` varchar(100) DEFAULT NULL,
`genome` varchar(250) DEFAULT NULL,
`nucleotide` text,
`aminoacid` text,
PRIMARY KEY (`id`),
KEY `gene_genome` (`seq_id`,`genome`)
)
''')
dna_folder = '/Volumes/Macintosh HD 2/thiberio/aeromonas_genomes/DNA_ORFs/formated'
aa_folder  = 'genomes_4_clustering'

dna_files = listdir(dna_folder)
aa_files = listdir(aa_folder)

already_in_db = []
cursor.execute('SELECT DISTINCT(genome) FROM sequences_all;')
for result in cursor.fetchall():
    already_in_db.append(result[0])

for tmp in aa_files:
    if not tmp.endswith('.faa'):
        continue

    genome_name = tmp.replace('.faa', '')

    if genome_name in already_in_db:
        print '\t**%s already on db!' %genome_name
        continue

    sequences = {}
    print genome_name

    fna = open('%s/%s.fna' %(dna_folder, genome_name)).read()
    for seq_block in fna.split('>'):
        if seq_block == '':
            continue
        seq_block = seq_block.strip()
        seq_block = seq_block.split('\n')
        header = re.search('^(\S+)\s\[(.+?)\]', seq_block[0], re.M)
        header = header.group(1)
        seq = ''.join(seq_block[1:]).upper()

        sequences[header] = {}
        sequences[header]['dna'] = seq
    del(fna)

    faa = open('%s/%s.faa' %(aa_folder, genome_name)) .read()
    for seq_block in faa.split('>'):
        if seq_block == '':
            continue
        seq_block = seq_block.strip()
        seq_block = seq_block.split('\n')
        header = re.search('^(\S+)\s\[(.+?)\]', seq_block[0], re.M)
        header = header.group(1)
        seq = ''.join(seq_block[1:]).upper()

        if header not in sequences:
            sequences[header] = {}
        sequences[header]['aa'] = seq
    del(faa)

    no_aa  = 0
    no_dna = 0
    sql_command = []
    for seq_id in sequences:
        if 'dna' not in sequences[seq_id]:
            no_dna += 1
        if 'aa' not in sequences[seq_id]:
            no_aa  += 1

        cursor.execute('INSERT INTO sequences_all (seq_id , genome , aminoacid, nucleotide) VALUES ("%s" , "%s" , "%s" , "%s")' %(seq_id, genome_name, sequences[seq_id]['aa'], sequences[seq_id]['dna']))
                
    if no_aa:
        print '**%i sequences have no aminoacid sequences' %no_aa
    if no_dna:
        print '**%i sequences have no nucleotide sequences' %no_dna
    print ''

db.commit()
