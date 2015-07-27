__author__ = 'Thiberio'

from pickle import load
import MySQLdb as sql

db  = sql.connect(user='root', db='aeromonas')
cursor = db.cursor()

t3ss_groups = load(open('t3ss_groups-merged.pkl'))
groups = load(open('homologous_groups-merged.pkl'))

for group in t3ss_groups:
    sql_query = []
    for genome in groups[group]:
        for protein in groups[group][genome]:
            sql_query.append('(genome="%s" and seq_id="%s")' %(genome, protein))
    cursor.execute('select seq_id, genome, aminoacid from sequences where %s' %' or '.join(sql_query))

    out = open('interproscan-merged/fastas/%s.faa' %group.replace('&', '-'), 'wb')
    for query_result in cursor.fetchall():
        out.write('>%s|%s\n%s\n' %query_result)
    out.close()