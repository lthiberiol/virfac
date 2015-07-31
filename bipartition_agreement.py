from Bio import Phylo as phy
import pandas as pd
from os import listdir, chdir

#test

count = 0
reliable_branches = set()
respective_groups = {}
species_in_groups = {}
for tree_file in listdir('fasttree'):

    if not tree_file.endswith('.tre'):
        continue

    group = tree_file.replace('.tre', '')

    tree = phy.read('fasttree/%s' %tree_file, 'newick')

    species_in_groups[group] = set(tree.get_terminals())

    respective_groups[group] = []

    for branch in tree.get_nonterminals():
        if branch.confidence >= 0.9:
            leaves = set()
            for leaf in branch.get_terminals():
                leaves.add( leaf.name.split('|')[1] )
        
            reliable_branches.add(frozenset(leaves))
            respective_groups[group].append( frozenset(leaves) )

    count += 1
    if count == 5:
        break

compatibility_df = pd.DataFrame(index=reliable_branches, columns=respective_groups.keys(), data=0)
for group in respective_groups:
    for comparison_group in respective_groups:

        if group == comparison_group:
            continue
        
        for branch in respective_groups[branch]:
            for comparison_branch in respective_groups[comparison_group]:
            
                branch_union        = branch.union(       comparison_branch)
                branch_intersection = branch.intersection(comparison_branch)

                if branch_union == branch or branch_union == comparison_branch:
                    compatible = True
                elif not branch_intersection:
                    compatible = True









