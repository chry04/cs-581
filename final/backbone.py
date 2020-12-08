import sys
import os
import utils
import dendropy
import shutil

if __name__ == "__main__":
    tree_path = str(sys.argv[1])
    output = str(sys.argv[2])
    aln = str(sys.argv[3])
    query = str(sys.argv[4])

    # output path, ref, query

    aln_dict = utils.read_data(aln)
    ref_dict, q_dict = utils.seperate(aln_dict, query)
    tree = dendropy.Tree.get(path=tree_path, 
    			      schema="newick",
    			      rooting= 'force-unrooted')
    tree.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    			      
    base_bipart = tree.encode_bipartitions()
    namespace = tree.taxon_namespace
        
    for name, seq in q_dict.items():
        subtree = dendropy.Tree(tree)
        #print(name)
        subtree.prune_taxa_with_labels([name])
        
        subtree.collapse_basal_bifurcation(\
        	      set_as_unrooted_tree=True)

        tmp_tree = output+"/{}/backbone2.tree".format(name)
        #print(tmp_tree)
        subtree.write(path=tmp_tree, schema="newick",\
                      suppress_rooting=True)
    
   
