import re
from os import listdir

stop_line = '-------------------------------------------------------------'

yeah={}
for search in listdir('uniref90/'):
    if not search.endswith('.hmmout'):
        continue
    print ''
    print search
    entry = open('uniref90/%s' %search).read()
    flag = False
    yeah[search.split('.hmmout')[0]] = {}
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

        if re.search('(hypothetical|uncharacterised|unkown)', hit_desc):
#            print '\t%.2e' %evalue
#            print '\t'+hit_desc
            continue

        yeah[search.split('.hmmout')[0]] = {'evalue':evalue, 'min evalue':min_evalue, 'his desc':hit_desc}
        break

no_aero={}
for search in listdir('uniref90/'):
    if not search.endswith('.hmmout'):
        continue
    print ''
    print search
    entry = open('uniref90/%s' %search).read()
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
