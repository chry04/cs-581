import sys
import os
import utils
import dendropy
import shutil

if __name__ == "__main__":
    model = str(sys.argv[1])
    info = str(sys.argv[2])
    tree_path = str(sys.argv[3])
    output = str(sys.argv[4])
    aln = str(sys.argv[5])
    query = str(sys.argv[6])
    n = int(sys.argv[7])

    # output path, ref, query, backbone tree, info

    aln_dict = utils.read_data(aln)
    ref_dict, q_dict = utils.seperate(aln_dict, query)
    tree = dendropy.Tree.get(path=tree_path, schema="newick")

    files = []

    os.mkdir("tmp")

    for name, seq in q_dict.items():
        y = utils.find_y(seq, ref_dict)
        nodes = utils.subtree_nodes(tree, y, n)
        subtree = tree.extract_tree_with_taxa(nodes)

        tmp_tree = "tmp/tree_"+name
        tmp_aln = "tmp/aln"+name+".fa"
        tmp_output = "tmp/"+name+".jplace"
        subtree.write(path=tmp_tree, schema="newick")

        f = open(tmp_aln, "w")
        f.write(">"+name)
        f.write("\n")
        f.write(seq+"\n")
        for n in nodes:
            f.write(">"+n.label+"\n")
            f.write(ref_dict[n.label])
            f.write("\n")

        f.close()

        os.system("./pplacer -m {} -s {} -t {} -o {} {}".format(model, info, tmp_tree, tmp_output,\
                tmp_aln))
        files.append(tmp_output)

    jplace = ""
    for file in files:
        jplace += " "
        jplace += file

    os.system("./guppy merge -o {}{}".format("out.jplace", jplace))

    #shutil.rmtree("tmp")
