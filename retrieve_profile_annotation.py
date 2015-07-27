__author__ = 'Thiberio'
import re
from os import listdir
from pickle import load

stop_line = '-------------------------------------------------------------'
target_folder = 'profile_annotation-merged/uniref90'
group_desc = load(open('homologous_groups_descriptions-merged.pkl'))

yeah={}
for search in listdir(target_folder):

    if not search.endswith('.hmmout'):
        continue

    print search+'\n'

    entry = open('%s/%s' %(target_folder, search)).read()
    flag = False
#    yeah[search.split('.hmmout')[0]] = {}
    min_evalue = 10
    for line in entry.split('\n'):

        if not flag and line.replace(' ', '') == stop_line:
            flag = True
            continue

        if not flag:
            continue

        line = line.split()
        evalue = float(line[0])
        hit_id = line[8]

        if evalue > 10e-5:
            break

        if evalue < min_evalue:
            min_evalue = evalue
            best_hit_tax = re.search('^>>\s%s.*?Tax=(.*?)$' %hit_id, entry, re.M).group(1)
            if 'RepID' in best_hit_tax:
                best_hit_tax = best_hit_tax.split('RepID')[0].strip()

        hit_desc = re.search('^>>\s%s\s*(.*?)\s+n=' %re.escape(hit_id), entry, re.M)
        hit_desc = hit_desc.group(1)
        hit_tax  = re.search('^>>\s%s.*?Tax=(.*?)$' %hit_id, entry, re.M).group(1)
        if 'RepID' in hit_tax:
                hit_tax = hit_tax.split('RepID')[0].strip()

        if re.search('(hypothetical|uncharacterized|unkown)', hit_desc, re.I):
#            print '\t%.2e' %evalue
#            print '\t'+hit_desc
            continue

        yeah[search.split('.hmmout')[0]] = {'evalue':evalue, 'min evalue':min_evalue, 'hit desc':hit_desc, 'hit tax': hit_tax, 'best hit tax': best_hit_tax}
        break

no_aero={}
for search in listdir(target_folder):
    if not search.endswith('.hmmout'):
        continue
    print ''
    print search
    entry = open('%s/%s' %(target_folder, search)).read()
    flag = False
    no_aero[search.split('.hmmout')[0]] = {}
    min_evalue = 10
    for line in entry.split('\n'):
        if not flag and line.replace(' ', '') == stop_line:
            flag = True
            continue
        if not flag:
            continue

        line = line.split()
        evalue = float(line[0])
        hit_id = line[8]

        if evalue > 10e-5:
            break

        if evalue < min_evalue:
            min_evalue = evalue

        hit_desc = re.search('^>>\s%s\s*(.*?)$' %re.escape(hit_id), entry, re.M)
        hit_desc = hit_desc.group(1)
        if 'gi|' in hit_desc:
            hit_desc = hit_desc.split('gi|')[0]

        if re.search('(hypothetical|uncharacterised|unkown|aeromonas)', hit_desc, re.I):
#            print '\t%.2e' %evalue
#            print '\t'+hit_desc
            continue

        no_aero[search.split('.hmmout')[0]] = {'evalue':evalue, 'min evalue':min_evalue, 'his desc':hit_desc}
        break

# generate HTML
out = open('profile_annotations-merged.html', 'wb')
out.write('''
<html>
<body>
<table border="1">
    <tr><th>group id</th>
    <th>VFDB's description</th>
    <th>best hit e-value</th>
    <th>best hit source</th>
    <th>e-value</th>
    <th>description</th>
    <th>source</th></tr>
''')
for group in yeah:
    out.write('<tr style="text-align:center"> <th style="max-width:15em; text-align:left; ">%s</th> <td>%s</td> <td>%.2e</td> <td style="max-width:10em;">%s</td> <td>%.2e</td> <td>%s</td> <td style="max-width:10em;">%s</td> </tr>\n'
              %(group,
                group_desc[group.replace('-', '&')],
                yeah[group]['min evalue'],
                yeah[group]['best hit tax'],
                yeah[group]['evalue'],
                yeah[group]['hit desc'],
                yeah[group]['hit tax']))
out.write('''
</table>
</body>
</html>
''')
out.close()