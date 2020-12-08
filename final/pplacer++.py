import sys
import os
import utils
import dendropy
import shutil
from dendropy.calculate import treecompare

if __name__ == "__main__":
    model = str(sys.argv[1])
    info = str(sys.argv[2])
    tree_path = str(sys.argv[3])
    output = str(sys.argv[4])
    aln = str(sys.argv[5])
    query = str(sys.argv[6])
    n = int(sys.argv[7])
    run = int(sys.argv[8])

    # output path, ref, query, backbone tree, info

    aln_dict = utils.read_data(aln)
    ref_dict, q_dict = utils.seperate(aln_dict, query)
    tree = dendropy.Tree.get(path=tree_path, schema="newick", rooting='force-unrooted')
    #base_bipart = tree.encode_bipartitions()
    #tree.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    #tree.update_bipartitions()
    namespace = tree.taxon_namespace

    files = []
    try: 
        os.mkdir("tmp{}".format(run))
    except:
    	print("tmp directory already exists")    
    try:
        os.mkdir(output)
    except:
    	print("path directory already exists")
        

    for name, seq in q_dict.items():
        y = utils.find_y(seq, ref_dict)
        nodes = utils.subtree_nodes(tree, y, n)
        subtree = tree.extract_tree_with_taxa(nodes)
        #subtree = tree.extract_tree_with_taxa(nodes, suppress_unifurcations=False)
        
        #subtree.update_bipartitions()
        #subtree.update_taxon_namespace()
        
        subtree.collapse_basal_bifurcation(set_as_unrooted_tree=True)
        #print(dendropy.calculate.treecompare.false_positives_and_negatives(tree,subtree))
        #print(treecompare.weighted_robinson_foulds_distance(tree,subtree))

        tmp_tree = "tmp{}/tree_".format(run)+name
        tmp_tree1 = "tmp{}/tree1_".format(run)+name
        tmp_aln = "tmp{}/aln".format(run)+name+".fa"
        tmp_output = "tmp{}/".format(run)+name+".jplace"
        ful_output = output+"/"+name+".jplace" 
        subtree.write(path=tmp_tree, schema="newick", suppress_rooting=True)
        
       

        f = open(tmp_aln, "w")
        f.write(">"+name)
        f.write("\n")
        f.write(seq+"\n")
        for n in nodes:
            f.write(">"+n.label+"\n")
            f.write(ref_dict[n.label])
            f.write("\n")

        f.close()

        os.system("./pplacer -m {} -s {} -t {} --keep-at-most 1 -o {} {}".format(model, info, tmp_tree, ful_output,\
                tmp_aln))

        
        os.system("./guppy tog -o {} {}".format(tmp_tree1, ful_output))
        
        query_taxon = dendropy.Taxon(name)
        namespace.add_taxon(query_taxon)

        added_tree = dendropy.Tree.get(path=tmp_tree1, 
        				schema="newick",
        				taxon_namespace=namespace)

        place_leaf = added_tree.find_node_with_taxon_label(name)
        internal = place_leaf.adjacent_nodes()[0]
        direction = []
        dir_edge = []
        if internal.parent_node!=place_leaf:
            direction.append(internal.parent_node)
            dir_edge.append(internal.edge_length)

        for child in internal.child_nodes():
            if child != place_leaf:
                direction.append(child)
                dir_edge.append(child.edge_length)

        left, path_l = utils.find_closest(direction[0], {place_leaf, internal, direction[0]})
        right, path_r = utils.find_closest(direction[1], {place_leaf, internal, direction[1]})

        left = tree.find_node_with_taxon_label(left.taxon.label)
        right = tree.find_node_with_taxon_label(right.taxon.label)
        _, path = utils.find_closest(left, {left}, y=right)

        length = sum([x.length for x in path_l])+dir_edge[0]
        target_edge = path[-1]

        for i in range(len(path)-1):
            length -= path[i].length
            if length < 0:
                target_edge = path[i]
                break


        if target_edge:
            tree.is_rooted = True
            tree.reroot_at_edge(target_edge, length1=-length, length2=target_edge.length+length,\
                    update_bipartitions=True, suppress_unifurcations=False)

            q_taxon = dendropy.Taxon(name)
            tree.seed_node.new_child(taxon=q_taxon, edge_length=place_leaf.edge_length)

            tree.taxon_namespace.add_taxon(q_taxon)
            tree.collapse_basal_bifurcation(set_as_unrooted_tree=True)


        
        #placement = tree.find_node_with_taxon_label(y)
        #parent_node = 

        
        
        '''
        added_bipart = added_tree.encode_bipartitions()
        tree.update_bipartitions()
        base_bipart = tree.encode_bipartitions()
        
        full_bipart = []
        full_bipart = added_bipart + base_bipart
        

        result_tree = dendropy.Tree.from_bipartition_encoding(full_bipart, namespace)
        '''


                
        #result_tree.collapse_basal_bifurcation(set_as_unrooted_tree=True)        
        #result_tree.suppress_unifurcations(update_bipartitions=True)
        #added_bipart = dendropy.calculate.treecompare.find_missing_bipartitions(added_tree, tree)        
        
        #result_tree = dendropy.Tree.from_bipartition_encoding(\
        #        base_bipart + added_bipart, namespace)
        tree.write(path = output+"/"+name+".tree", schema="newick", suppress_rooting=True) 
        #namespace.remove_taxon(query_taxon)
          
    
    #shutil.rmtree("tmp{}".format(run))
