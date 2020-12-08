from dendropy import *
import heapq

#separete the query and ref sequence from the alignment file

def read_data(aln):
    f = open(aln)
    result = dict()

    taxa = ""
    seq = ""
    for line in f:
        if line[0] == '>':
            if taxa != "":
                result[taxa] = seq
            taxa = line[1:-1]
            seq = ""

        elif line == "/n":
            continue
        else:
            seq += line[:-1]
            
    if taxa != "":
        result[taxa] = seq


    return result

def seperate(aln_dict, query):
    f = open(query)
    q_name = set()
    for line in f:
        q_name.add(line[:-1])

    print(q_name)

    ref = dict()
    query = dict()

    for key, value in aln_dict.items():
        if key in q_name:
            query[key] = value
        else:
            ref[key] = value
    
    return ref, query

def hamming(seq1, seq2):
    return len([1 for i in range(len(seq1)) if seq1[i] != seq2[i]])


def find_y(x, ref):
    low = len(x)
    y = ""
    for name, seq in ref.items():
        h_dist = hamming(x, seq)
        if h_dist < low:
            low = h_dist
            y = name
    return y


def subtree_nodes(tree, y, n):
    leaf_y = tree.find_node_with_taxon_label(y)
    #print(leaf_y)
    queue = [(0, 0, leaf_y.parent_node)]
    #print(queue)
    leaves = []
    visited = {y}

    counter = 1

    while len(leaves) < n:
        try:
            (length, _, node) = heapq.heappop(queue)
        except IndexError:
            break

        visited.add(node)
        if node.is_leaf():
            leaves.append(node)

        for child in node.adjacent_nodes():
            if child not in visited:
                heapq.heappush(queue, (length+1, counter, child))
                counter += 1

    
    result = []
    for item in leaves:
        result.append(item.taxon)

    return result


def compareTreesFromPath(treePath1, treePath2):
    print("Comparing {} with {}".format(treePath1, treePath2))

    tax = TaxonNamespace()
    tr1 = Tree.get(path=treePath1,
                            schema='newick',
                            rooting='force-unrooted',
                            taxon_namespace=tax,
                            preserve_underscores=True)
    tr2 = Tree.get(path=treePath2,
                            schema='newick',
                            rooting='force-unrooted',
                            taxon_namespace=tax,
                            preserve_underscores=True)

    tr1.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    tr2.collapse_basal_bifurcation(set_as_unrooted_tree=True)

    return compareDendropyTrees(tr1, tr2)
    # print("RF distance on %d shared leaves: %d" % (nl, fp + fn))


def compareDendropyTrees(tr1, tr2):
    from dendropy.calculate.treecompare \
        import false_positives_and_negatives

    lb1 = set([l.taxon.label for l in tr1.leaf_nodes()])
    lb2 = set([l.taxon.label for l in tr2.leaf_nodes()])

    com = lb1.intersection(lb2)
    if com != lb1 or com != lb2:
        com = list(com)
        tns = TaxonNamespace(com)

        tr1.retain_taxa_with_labels(com)
        tr1.migrate_taxon_namespace(tns)

        tr2.retain_taxa_with_labels(com)
        tr2.migrate_taxon_namespace(tns)
    com = list(com)

    tr1.update_bipartitions()
    tr2.update_bipartitions()

    nl = len(com)
    ei1 = len(tr1.internal_edges(exclude_seed_edge=True))
    ei2 = len(tr2.internal_edges(exclude_seed_edge=True))

    [fp, fn] = false_positives_and_negatives(tr1, tr2)
    rf = float(fp + fn) / (ei1 + ei2)

    return (nl, ei1, ei2, fp, fn, rf)
