import ete3
from Bio import Phylo as phy
import os
import re

class cd: 
    """
    Context manager for changing the current working directory
    """
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

reconciliations = open('ranger.output_final-gene_names').read()
group_order     = re.findall( '^ ------------ Reconciliation for Gene Tree (\S+) \(Unrooted\) -------------$', reconciliations, re.M )
reconciliations = re.split( '^ ------------ Reconciliation for Gene Tree \S+ \(Unrooted\) -------------$', reconciliations, flags=re.M )
reconciliations.pop(0)

def root_like( ref_tree, tree_to_root ):
    outgroup = ''
    for node in sorted( ref_tree.children, key=len ):
        if node.is_leaf():
            putative_outgroups = tree_to_root.get_leaves_by_name(node.name)
            for leaf in putative_outgroups:
                tree_to_root.set_outgroup( leaf )
                if tree_to_root.get_topology_id() == ref_tree.get_topology_id():
                    outgroup = leaf
                    break
        else:
            outgroup_members = []
            for leaf in node:
                outgroup_members.append( tree_to_root.get_leaves_by_name(leaf.name) )

            for outgroup_combination in itertools.product( *outgroup_members ):
                if tree_to_root.check_monophyly( [x.full_name for x in outgroup_combination], 'full_name' )[0]:
                    putative_outgroup = tree_to_root.get_common_ancestor( outgroup_combination )
                    tree_to_root.set_outgroup( putative_outgroup )
                    if tree_to_root.get_topology_id() == ref_tree.get_topology_id():
                        outgroup = putative_outgroup
                        break
        if outgroup:
            return tree_to_root

    if not outgroup:
        return None

def name_matching_branches( ref_tree, unamed_tree ):
    for node in ref_tree.traverse():
        if node.is_leaf():
            continue

        branch_members = []
        for leaf in node:
            branch_members.append( unamed_tree.get_leaves_by_name(leaf.name) )

        for members_combination in itertools.product( *branch_members ):
            putative_match = unamed_tree.get_common_ancestor( members_combination )
            if node.get_topology_id() == putative_match.get_topology_id() and set(node.get_leaf_names()) == set(putative_match.get_leaf_names()):
                putative_match.name = node.name
                break
    return unamed_tree.copy()


without_branch_names = []
with_branch_names    = []
times                = []
for group, reconciliation in zip( group_order, reconciliations ):
    print group
    recon_tree = ete3.Tree( reconciliation.split('\n')[5], format=1)
    original_tree = ete3.Tree( '%s/RAxML_bipartitions.%s-final' %(group, group) )
    for node in original_tree.traverse():
        if node.is_leaf() and '|' in node.name:
            if node.name.startswith('A_'):
                node.full_name = node.name
                node.name      = node.name.split('|')[0]
            else:
                node.full_name = node.name
                node.name      = node.name.split('|')[1]
        elif node.is_leaf() and '|' not in node.name:
            node.full_name = node.name

    start_time = time()
    rooted_tree = root_like( recon_tree, original_tree )
    if rooted_tree:
        final = name_matching_branches( recon_tree, rooted_tree )
        times.append( time() - start_time )
        for i in final.traverse():
            if i.is_leaf():
                continue

            j = recon_tree.search_nodes( name=i.name )[0]
            if j.get_topology_id() != i.get_topology_id():
                print n.name

    else:
        print "\tFail"

print 'yeah'
