import re
from os.path import getsize, isfile
from Bio.SearchIO import parse
from pickle import load, dump
from os import system, listdir, chdir
import multiprocessing

groups = load( open( 'homologous_groups-merged.pkl' ) )
t6ss_groups = load( open( 't6ss_groups-merged.pkl' ) )

system('mkdir profile_annotation-merged/fastas')
for group in t6ss_groups:
    if len(group) > 200:
        group_name = group[:200]
    else:
        group_name = group
    if '-' in group:
        system('cat ../homologous_groups/%s.faa > profile_annotation-merged/fastas/%s.faa' %('.faa ../homologous_groups/'.join(group.split('-') ), group_name ) )
    else:
        system('cp ../homologous_groups/%s.faa profile_annotation-merged/fastas/' %group)

system('mkdir profile_annotation-merged/mafft')
def run_mafft(fasta):
    system('mafft --reorder --auto profile_annotation-merged/fastas/%s > profile_annotation-merged/mafft/%s' %( fasta, fasta.replace( '.faa', '.aln' ) ) )
pool = multiprocessing.Pool(processes = 10)
pool.map(run_mafft, listdir('profile_annotation-merged/fastas'))
pool.close()
pool.join()

system('mkdir profile_annotation-merged/hmm')
def run_hmmbuild(aln):
    if getsize('profile_annotation-merged/mafft/%s' %aln):
       system('hmmbuild profile_annotation-merged/hmm/%s profile_annotation-merged/mafft/%s' %(aln.replace('.aln', '.hmm'), aln) ) 
    else:
       system('hmmbuild profile_annotation-merged/hmm/%s profile_annotation-merged/fastas/%s' %(aln.replace('.aln', '.hmm'), aln.replace('.aln', '.faa')) ) 
pool = multiprocessing.Pool(processes = 10)
pool.map(run_hmmbuild, listdir('profile_annotation-merged/mafft'))
pool.close()
pool.join()
    
system('mkdir profile_annotation-merged/hmmout')
def run_hmmsearch(profile):
    system('hmmsearch --cpu 2 -o profile_annotation-merged/hmmout/%s profile_annotation-merged/hmm/%s ~/work/blast_dbs/uniref90.fasta' %( profile+'out', profile) )
pool = multiprocessing.Pool(processes = 10)
pool.map(run_hmmsearch, listdir('profile_annotation-merged/hmm'))
pool.close()
pool.join()


stop_line = '-------------------------------------------------------------'
target_folder = 'profile_annotation-merged/hmmout'
group_desc = load(open('homologous_groups_descriptions-merged.pkl'))

group_annotation={}
for hmmout in listdir(target_folder):

    if not hmmout.endswith('.hmmout'):
        continue

    print hmmout+'\n'

    entry = open('%s/%s' %(target_folder, hmmout)).read()
    flag = False
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
            continue

        group_annotation[hmmout.split('.hmmout')[0]] = {'evalue':evalue, 'min evalue':min_evalue, 'hit desc':hit_desc, 'hit tax': hit_tax, 'best hit tax': best_hit_tax}
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
for group in group_annotation:

    for group_name in group_desc:
        if group in group_name:
            break

    out.write('<tr style="text-align:center"> <th style="max-width:15em; text-align:left; ">%s</th> <td>%s</td> <td>%.2e</td> <td style="max-width:10em;">%s</td> <td>%.2e</td> <td>%s</td> <td style="max-width:10em;">%s</td> </tr>\n'
              %(group,
                group_desc[group_name],
                group_annotation[group]['min evalue'],
                group_annotation[group]['best hit tax'],
                group_annotation[group]['evalue'],
                group_annotation[group]['hit desc'],
                group_annotation[group]['hit tax']))
out.write('''
</table>
</body>
</html>
''')
out.close()
