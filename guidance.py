#!/usr/bin/python
# -*- coding: utf-8 -*-

import pylab as pl
import pandas as pd 
from pickle import load
from sys import exit
from os import system
import re
import multiprocessing
from multiprocessing import Manager
manager = Manager()
num_threads = 20

hg = open('homologous_groups.list').readlines()
target_clusters = []
min_size = 2

def run_guidance(start, end):
    for element in target_clusters[start:end]:
        print element
        system('nice mafft --auto homologous_groups/%s.faa > aligned_groups_all2/%s.aln' %(element, element))

for line in hg:
    if line.startswith(':'):
        continue
    regex = re.search('^cluster\s(.*?)\ssize=(\d+?)\s', line, re.M)
    cluster_name = regex.group(1)
    cluster_size = int(regex.group(2))
    if cluster_size < min_size:
        continue
    target_clusters.append(cluster_name)

jobs = []
range_threads = len(target_clusters)/num_threads
for thread in range(0, num_threads):
    start = thread*range_threads
    end = start + range_threads
    print '**Creating thread %i!' %thread
    process = multiprocessing.Process(target=run_guidance, args=(start, end))
    process.start()
    jobs.append(process)
for p in jobs:
    p.join()
