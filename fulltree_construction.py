#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run forest output, needs IQTREE2
"""

from forest_functions import *

## TRUE distance matrix example test to reconstruct full tree by varying input parameters
#np.random.seed(2024)
iqtree2_path = "./qtree2"

Ntips = [32]
lb = 0.05
ub = 0.1
Thr = 2
M_steps = 20
m_steps = 40
taus = np.linspace(0.01, lb/2+0.02, 30) #sample more closer to 0 less to 0.1
#taus = np.linspace(lb/2, ub, 100) #sample more closer to 0 less to 0.1
#exponents = np.linspace(-5, 0, 50)  
#taus = 0.1 * np.exp(exponents)
print(f"taus: {taus}")
epsilon = 1e-8
Nruns = 5
#replicat = 1
k_values = [1024]#[64, 128, 256, 512, 1024, 2048, 4096]
num_aln_reps = 2 #should be greater than 1 or otherwise IQTREE will not produce correct file name
    
for ntips in Ntips:
    tipdir = f"./rand_nj_lowestrate/{ntips}tips"
    shutil.rmtree(tipdir)
    os.makedirs(tipdir)
    print(tipdir)
    true_tree = dendropy.simulate.treesim.birth_death_tree(birth_rate=1.0, death_rate=0, num_extant_tips=ntips)
    for edge in true_tree.edges():
        edge.length = np.random.uniform(lb, ub)
    print(true_tree)
    for idx, leaf in enumerate(true_tree.leaf_nodes()):
        leaf.taxon.label = str(idx)
#     file_path = "./data/8tips/" + tree_shape + '_' + str(ntips) + ".tre"
#     true_tree = Phylo.read(file_path, 'newick')
#     for i, clade in enumerate(true_tree.find_clades(terminal=True)):
#         clade.name = str(i)
#     print("True tree")
#     Phylo.draw_ascii(true_tree)
#     Phylo.draw(true_tree, branch_labels = lambda c: c.branch_length)
#     true_newick_tree = io.StringIO()
#     Phylo.write(true_tree, true_newick_tree, "newick")
#     true_newick_string = true_newick_tree.getvalue()
#     print(f'True tree: {true_newick_tree.getvalue()}')
#     taxon_namespace = dendropy.TaxonNamespace()
#     true_tree = dendropy.Tree.get(data=true_newick_tree.getvalue(), schema="newick", taxon_namespace=taxon_namespace)

    tree_path = os.path.join(tipdir, f"true_tree_{ntips}.tre")
    with open(tree_path, 'w') as tree_file:
       tree_file.write(true_tree.as_string(schema="newick"))
    draw_true_tree = Phylo.read(tree_path, 'newick')
    Phylo.draw(draw_true_tree, branch_labels = lambda c: np.round(c.branch_length, 4))
    taxon_namespace = dendropy.TaxonNamespace()
    true_tree_string = re.search(r' (?=\()(.+;)', true_tree.as_string(schema="newick")).group(1) + '\n'
    true_tree_string = adjust_taxon_labels(true_tree_string)
    true_tree = dendropy.Tree.get(data=true_tree_string, schema="newick", taxon_namespace=taxon_namespace)
    pdm = dendropy.calculate.treemeasure.PatristicDistanceMatrix(true_tree)
    labels = [taxon.label for taxon in pdm.taxon_namespace]
    matrix_size = len(pdm.taxon_namespace)
    distance_matrix = [[0 for _ in range(matrix_size)] for _ in range(matrix_size)]
    for i, taxon1 in enumerate(pdm.taxon_namespace):
        for j, taxon2 in enumerate(pdm.taxon_namespace):
            distance_matrix[i][j] = pdm(taxon1, taxon2)
    true_distance_matrix = pd.DataFrame(distance_matrix, index=labels, columns=labels)
    true_distance_matrix.index = adjust_labels(true_distance_matrix.index)
    true_distance_matrix.columns = adjust_labels(true_distance_matrix.columns)
    true_distance_matrix = true_distance_matrix.sort_index().sort_index(axis=1)
    true_distance_matrix.index = true_distance_matrix.index.astype(int)  # Convert index to integer
    true_distance_matrix = true_distance_matrix.sort_index()  # Sort DataFrame by index
    true_distance_matrix.index = true_distance_matrix.index.astype(str)  

    true_distance_matrix.columns = true_distance_matrix.columns.astype(int)  # Convert column names to integer
    true_distance_matrix = true_distance_matrix.sort_index(axis=1)  # Sort DataFrame by column names
    true_distance_matrix.columns = true_distance_matrix.columns.astype(str)  # Convert column names back to string if necessary

    
for ntips in Ntips:
    a = true_distance_matrix
    print("Input Distance Matrix with a chord depth of ", calculate_chord_depth(true_tree_string, a))
    M = np.ceil(calculate_chord_depth(true_tree_string, a)) + 1
    internal_edges = get_internal_edge_lengths(Tree(true_tree_string))
    print(f"internal_edges: {0.5*np.sort(internal_edges)}")
    print("M, max tau: ", M, 0.5*np.sort(internal_edges)[1])
    epsilon = 1e-5
    Thr = ntips
#     epsilon = 1e-12
#     tau_step = 10
#     tau_values = np.linspace(epsilon, 0.5*np.sort(internal_edges)[1] + 0.01, tau_step)
#     tau_values = np.arange(epsilon, 0.5*np.sort(internal_edges)[1] + 0.01, 0.001)
#     combined_tau_list = np.union1d(tau_values, np.array([0.5*np.sort(internal_edges)[1] - epsilon, 0.5*np.sort(internal_edges)[1], 0.5*np.sort(internal_edges)[1]]))
#     tau_list = list(combined_tau_list)
    tau_list = list(taus)
    print("len(tau_list): ", len(tau_list))
    print("tau_list: ", tau_list)



    true_g = construct_clustering_graph(np.array(a), M)
    components = list(nx.connected_components(true_g))
    print("Number of components:", len(components))
    for i, component in enumerate(components):
        print("Component", i+1, "has", len(component), "nodes")
    random_subgraph = get_random_subgraph(true_g, Thr)


    tree_shape = 'bal'

    rf_dist_list = []
    leaves_pairs = [(u, v) for u in random_subgraph.nodes() for v in random_subgraph.nodes() if u < v]
    for tau in tau_list:
        bipartitions = []
        final_bipartitions = []
        mini_bipartitions_list = []
        for leaves in leaves_pairs:
            mini_bipartitions, clusters = mini_contractor(random_subgraph, np.array(a),  leaves, M, tau)
            mini_bipartitions_list.append(mini_bipartitions)
        flattened_bipartitions = [tuple(frozenset(part) for part in bipartition) for sublist in mini_bipartitions_list for bipartition in sublist]
        unique_bipartitions = list(set(flattened_bipartitions))
        unique_bipartitions_as_sets = [tuple(set(part) for part in bipartition) for bipartition in unique_bipartitions]
        unique_bipartitions_as_sets.sort()
        #print(f'{len(flattened_bipartitions)} bipartitions')
        sorted_frozenset_data = [sort_and_freeze(pair) for pair in unique_bipartitions_as_sets]
        unique_sorted_frozenset_data = set(sorted_frozenset_data)
        unique_sorted_tuples = [tuple(map(set, pair)) for pair in unique_sorted_frozenset_data]
        #print(f'{len(unique_sorted_tuples)} unique bipartitions')
        unique_trees = get_unique_trees(unique_sorted_tuples, ntips, 10)
        if len(unique_trees) == 1: #look for a single reconstructed tree
            for tree in unique_trees:
                print(f'tau: {tau}, 2*tau: {2*tau}, 4*tau: {4*tau}')
                print(tree)
                tree = Phylo.read(io.StringIO(tree), "newick")
                pd_tree = io.StringIO()
                Phylo.write(tree, pd_tree, "newick")
                pd_tree = re.sub(r'Inner\d+:', '', pd_tree.getvalue())
                pd_tree = dendropy.Tree.get(data=pd_tree, schema="newick", taxon_namespace=taxon_namespace)
                rf_dist = dendropy.calculate.treecompare.symmetric_difference(true_tree, pd_tree)
                print("RF distance: ", rf_dist)
                rf_dist_list.append(rf_dist)
        else:
            rf_dist_list.append(1000)

    print("tau: ", tau_list)
    print("rf: ", rf_dist_list)
    
tree_shape = "rand"
index_last_zero = max(loc for loc, val in enumerate(rf_dist_list) if val == 0)
unique_vals = np.unique(rf_dist_list)
second_largest = unique_vals[-2] if len(unique_vals) > 1 else unique_vals[0]
new_value = second_largest + 2
rf_dist_list = [new_value if x == 1000 else x for x in rf_dist_list]
plt.plot(tau_list, rf_dist_list, marker='o')
plt.axvline(x=np.sort(internal_edges)[0] * 0.5, color='red', linestyle='-', label=f'0.5 Min Internal Edge = {np.round(np.sort(internal_edges)[0] * 0.5, 3)}')
plt.axvline(x=tau_list[index_last_zero], color='green', linestyle='-', label=f'Largest Tau where RF is 0: {np.round(tau_list[index_last_zero], 3)}')
plt.yticks(list(np.unique(rf_dist_list)), [f'{int(x)}' if x != new_value else 'NA' for x in np.unique(rf_dist_list)])
plt.xlabel('Tau Values')
plt.ylabel('RF Distance')
plt.title(f'RF Distance vs Tau Values (True dist matrix, rand_{ntips} tips)')
plt.legend()
plt.savefig("./rand32_true.png")
plt.show()

print(f"1/2 internal lengths: {0.5 * np.sort(internal_edges)}")
#(0,3,(1,2),(6,27,28,29,30,31,(7,11,12,15,(13,14),(8,9,10)),(16,17,21,22,23,24,(20,(18,19))),(25,26),(4,5)));


#Est dist mat test
for k in k_values:
    print(f"seq_length: {k}")
    alnfile = os.path.join(tipdir, f"aln_{k}")
    print(alnfile, tree_path)
    cmd = [iqtree2_path, "--alisim", alnfile, "-t", tree_path, "-m", "JC", "--length", str(k), "--num-alignments", str(num_aln_reps), "--out-format", "fasta", "--seed", "2024"]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    for i_replicate in range(1, num_aln_reps+1):
        print(f"replicate: {i_replicate}")
        alignment_path = os.path.join(tipdir, f"aln_{k}_{i_replicate}.fa")
        alignments = SeqIO.parse(alignment_path, 'fasta')

        msa = MultipleSeqAlignment(alignments)
        calculator = DistanceCalculator('identity')
        est_dist_mat_orig = calculate_jc69_distance_matrix(msa)

        
        matrix_data = est_dist_mat_orig.matrix
        labels = est_dist_mat_orig.names
        matrix_data = est_dist_mat_orig.matrix
        full_matrix = [[0] * len(labels) for _ in range(len(labels))]
        for i in range(len(labels)):
            for j in range(i + 1):
                full_matrix[i][j] = matrix_data[i][j]
                full_matrix[j][i] = matrix_data[i][j]  

        est_dist_mat = pd.DataFrame(full_matrix, index=labels, columns=labels)


        constructor = DistanceTreeConstructor()
        tree = constructor.nj(est_dist_mat_orig)
        nj_tree = io.StringIO()
        Phylo.write(tree, nj_tree, "newick") 
        nj_tree_string = re.sub(r'Inner\d+:', '', nj_tree.getvalue())
        nj_tree_string = adjust_taxon_labels(nj_tree_string)
        nj_tree = dendropy.Tree.get(data=nj_tree_string.rstrip('\n'), schema="newick", taxon_namespace=taxon_namespace)
        rf_distance = dendropy.calculate.treecompare.symmetric_difference(true_tree, nj_tree)
        if rf_distance != 0:
            print(true_tree_string)
            print(nj_tree_string)
            print("Max distances difference d - d-hat")
            print(np.abs(true_distance_matrix - est_dist_mat).max().max())
#             print(f"ntips: {ntips}, replicat: {i_replicate}, rf: {rf_distance}")
#             est_dist_mat_file = os.path.join(tipdir, f"est_dist_mat_{k}_{i_replicate}.csv")
#             forest_results  = get_forest_m(true_tree_string, nj_tree_string, true_distance_matrix, est_dist_mat, m_steps, taus, M_steps, ntips, Thr, k, epsilon, Nruns, printout = 0)
#             forest_results


true_dist = true_distance_matrix
internal_edges = get_internal_edge_lengths(Tree(true_tree_string))
print(f"internal_edges: {np.sort(internal_edges)}")

print("M, max tau: ", M, 0.5*np.sort(internal_edges)[1])
epsilon = 1e-5
Thr = ntips
tau_list = list(taus)
print("len(tau_list): ", len(tau_list))
print("tau_list: ", tau_list)
print(f'diff < tau ', np.abs(true_distance_matrix - est_dist_mat).max().max())
print(f'd < M + tau ', true_distance_matrix.max().max())
print(f'd-hat < M + tau ', est_dist_mat.max().max())


true_g = construct_clustering_graph(np.array(true_dist), M)
components = list(nx.connected_components(true_g))
print("Number of components:", len(components))
for i, component in enumerate(components):
    print("Component", i+1, "has", len(component), "nodes")
random_subgraph = get_random_subgraph(true_g, Thr)


rf_dist_list = []
leaves_pairs = [(u, v) for u in random_subgraph.nodes() for v in random_subgraph.nodes() if u < v]
for tau in tau_list:
    bipartitions = []
    final_bipartitions = []
    mini_bipartitions_list = []
    for leaves in leaves_pairs:
        mini_bipartitions, clusters = mini_contractor(random_subgraph, np.array(est_dist_mat),  leaves, M, tau)
        mini_bipartitions_list.append(mini_bipartitions)
    flattened_bipartitions = [tuple(frozenset(part) for part in bipartition) for sublist in mini_bipartitions_list for bipartition in sublist]
    unique_bipartitions = list(set(flattened_bipartitions))
    unique_bipartitions_as_sets = [tuple(set(part) for part in bipartition) for bipartition in unique_bipartitions]
    unique_bipartitions_as_sets.sort()
    #print(f'{len(flattened_bipartitions)} bipartitions')
    sorted_frozenset_data = [sort_and_freeze(pair) for pair in unique_bipartitions_as_sets]
    unique_sorted_frozenset_data = set(sorted_frozenset_data)
    unique_sorted_tuples = [tuple(map(set, pair)) for pair in unique_sorted_frozenset_data]
    #print(f'{len(unique_sorted_tuples)} unique bipartitions')
    unique_trees = get_unique_trees(unique_sorted_tuples, ntips, 10)
    if len(unique_trees) == 1: #look for a single reconstructed tree
        for tree in unique_trees:
            print(f'tau: {tau}, M: {M}')
            print(tree)
            tree = Phylo.read(io.StringIO(tree), "newick")
            pd_tree = io.StringIO()
            Phylo.write(tree, pd_tree, "newick")
            pd_tree = re.sub(r'Inner\d+:', '', pd_tree.getvalue())
            pd_tree = dendropy.Tree.get(data=pd_tree, schema="newick", taxon_namespace=taxon_namespace)
            rf_dist = dendropy.calculate.treecompare.symmetric_difference(true_tree, pd_tree)
            print("RF distance: ", rf_dist)
            rf_dist_list.append(rf_dist)
    else:
        rf_dist_list.append(1000)

print("tau: ", tau_list)
print("rf: ", rf_dist_list)


index_last_zero = max(loc for loc, val in enumerate(rf_dist_list) if val == 0)
unique_vals = np.unique(rf_dist_list)
second_largest = unique_vals[-2] if len(unique_vals) > 1 else unique_vals[0]
new_value = second_largest + 2
rf_dist_list = [new_value if x == 1000 else x for x in rf_dist_list]
plt.plot(tau_list, rf_dist_list, marker='o')
plt.axvline(x=np.sort(internal_edges)[0] * 0.5, color='red', linestyle='-', label=f'0.5 Min Internal Edge = {np.round(np.sort(internal_edges)[0] * 0.5, 3)}')
plt.axvline(x=tau_list[index_last_zero], color='green', linestyle='-', label=f'Largest Tau where RF is 0: {np.round(tau_list[index_last_zero], 3)}')
plt.yticks(list(np.unique(rf_dist_list)), [f'{int(x)}' if x != new_value else 'NA' for x in np.unique(rf_dist_list)])
plt.xlabel('Tau Values')
plt.ylabel('RF Distance')
plt.title(f'RF Distance vs Tau Values (Est dist matrix, rand_{ntips} tips)')
plt.legend()
plt.savefig("./rand32_est.png")
plt.show()
print(f"1/2 internal lengths: {0.5 * np.sort(internal_edges)}")
