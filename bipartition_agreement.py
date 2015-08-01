from Bio import Phylo as phy
import pandas as pd
from os import listdir, chdir
from copy import deepcopy

tree_folder = 'fasttree'
count = 0
reliable_branches = set()
respective_groups = {}
species_in_groups = {}
for tree_file in listdir(tree_folder):

    print tree_file

    if not tree_file.endswith('.tre'):
        continue

    group = tree_file.replace('.tre', '')

    tree = phy.read('%s/%s' %(tree_folder, tree_file), 'newick')

    species_in_groups[group] = set()

    respective_groups[group] = []

    for branch in tree.get_nonterminals():
        if branch.confidence >= 0.9:
            leaves = set()
            for leaf in branch.get_terminals():
                leaves.add( leaf.name.split('|')[1] )
 
            reliable_branches.add( frozenset( leaves ) )
            respective_groups[group].append( frozenset( leaves ) )
            species_in_groups[group].update( leaves )

    count += 1
    if count == 5:
        break

indexes = []
for branch in reliable_branches:
    indexes.append( '--'.join( branch ) )
compatibility_df = pd.DataFrame(index=indexes, columns=respective_groups.keys(), data=0)

#
# looping through groups, groups will be compared to all other branches
#   (not from the same group);
#
for group in respective_groups:

    #
    # looping through ALL BRANCHES!
    #
    for branch in reliable_branches:

        index = '--'.join( branch )
        
        #
        # if branch is part of the set of branches of a group, it is compatible
        #   so continue as nothing had happened
        #
        if branch in respective_groups[group]:
            continue

        #
        # test if there is species in the branch that are not present in the set
        #   of species from the group. If there is no intersection branch is
        #   compatible.
        #
        branch_species_not_in_group = branch.difference( species_in_groups[group] )
        if len(branch) == len(branch_species_not_in_group):
            continue

        #
        # if there are species on branch there are not present in group, remove
        #   them...
        #
        ref_branch = set( deepcopy( branch ) )
        for removable_species in branch_species_not_in_group:
            ref_branch.remove( removable_species )

        if len(ref_branch) < 2:
            continue

        #
        # compare the reference branch against all branches from group. In case
        #   of incompatibility, set respective dataframe position to 1...
        #
        for comparison_branch in respective_groups[group]:
            
            if ref_branch.isdisjoint(comparison_branch):
                continue
            if ref_branch.issubset(comparison_branch) or ref_branch.issuperset(comparison_branch):
                continue

            comparison_branch_complement = species_in_groups[group].difference(comparison_branch)
            if ref_branch.isdisjoint(comparison_branch_complement):
                continue
            if ref_branch.issubset(comparison_branch_complement) or ref_branch.issuperset(comparison_branch_complement):
                continue

            compatibility_df.loc[index][group] = 1
            break
