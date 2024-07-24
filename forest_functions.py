#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
functions to run forest results

"""

import numpy as np
import networkx as nx
import random
import dendropy
from dendropy import Tree, TaxonNamespace, Bipartition
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
import re
import io
from io import StringIO
import os
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
np.random.seed(0)


def construct_clustering_graph(distance_matrix, M):
    num_leaves = distance_matrix.shape[0]
    G = nx.Graph()
    for i in range(num_leaves):
        G.add_node(i)
    for i in range(num_leaves):
        for j in range(i+1, num_leaves):
            if distance_matrix[i, j] < M:
                G.add_edge(i, j)
    return G


def get_random_subgraph(graph, min_nodes):
    clusters = [graph.subgraph(c).copy() for c in nx.connected_components(graph)]
    suitable_clusters = [c for c in clusters if c.number_of_nodes() >= min_nodes]
    if not suitable_clusters:
        return None
    return random.choice(suitable_clusters)


def sort_and_freeze(pair):
    return tuple(sorted((frozenset(s) for s in pair), key=lambda s: (min(s) if s else float('inf'))))


def reconstruct_unrooted_from_bipartitions(unique_sorted_tuples, num_taxa):
    number_strings = [str(i) for i in range(num_taxa)]
    labels_names = dendropy.TaxonNamespace(number_strings, label="taxa1")
    split_matrix = np.zeros((len(unique_sorted_tuples), num_taxa))

    for i, (A, B) in enumerate(unique_sorted_tuples):
        for a in A:
            split_matrix[i, a] = 1

    #print(split_matrix.astype(int))
    bitmask_splits = []
    for bipartition in split_matrix:
        bitmask = 0
        for index, value in enumerate(bipartition):
            if value == 1:
                bitmask |= 1 << index
        bitmask_splits.append(bitmask)
    #print(bitmask_splits)
    number_strings = [str(i) for i in range(num_taxa)]
    my_unrooted_tree = dendropy.Tree.from_split_bitmasks(bitmask_splits, taxon_namespace = labels_names)
    #print(my_unrooted_tree.as_string(schema="newick"))
    #my_rooted_tree = dendropy.Tree.from_split_bitmasks(bitmask_splits, taxon_namespace = labels_names, is_rooted = True)
    #print(my_rooted_tree.as_string(schema="newick"))
    return my_unrooted_tree


#Forest algorithm functions
#m < 0.5(M - 3tau); 
#Typically, M is much larger than tau, so reconstructs a subforest of T with chord depth ≈1/2M which 
#all edges < tau will be contracted

# Phi is a distance from node u to the intersection point of (u, v, w)
def calculate_phi(u, v, w, distance_matrix):
    return 0.5 * (distance_matrix[u][v] + distance_matrix[u][w] - distance_matrix[v][w])

#works for pair of leaves in the given connected graph component, result depends on M, tau  
def mini_contractor(component, distance_matrix, leaves, M, tau):
    u, v = leaves
    B = {w for w in component if max(distance_matrix[u][w], distance_matrix[v][w]) < M}
    if u not in B:
        B.add(u)  
    if v not in B:
        B.add(v)  
    #print("For a pair of leaves: ", u, " and ", " v ", v, " the ball is: ", B) 
    S = B - {u}
    bipartitions = []
    x_minus_1 = u
    j = 0
    C = {0: [u]} 
    phi = {w: calculate_phi(u, v, w, distance_matrix) for w in S}
    #print("Phi from mini-contractor for each node w: ", phi)
    while S:
        x_0 = min(S, key=lambda w: phi[w])
        if j == 0: 
            phi[u] = 0
        #print("Phi(x_0), Phi(x_minus_1), 2*tau: ", phi[x_0], phi[x_minus_1], 2*tau)
        if phi[x_0] - phi[x_minus_1] >= 2 * tau:
        #print("long edge (bipartition) created between: ", x_0, " and ", x_minus_1)
            bipartitions.append([B - S.copy(), S.copy()])
            C[j + 1] = [x_0]
            j += 1
        else:
        #print(x_0, " is close to ", x_minus_1, " no bipartition is created")
            C[j].append(x_0)
        #print("Bipartitions and j: ", bipartitions, j)
        #print("Clusters: ", C)
        # Update the set S (reducing) and the previous node
        S.remove(x_0)
        x_minus_1 = x_0
    return bipartitions, {k: list(v) for k, v in C.items()}


# extends bipartitions in cases if there are other nodes in the graph component but not in bipartitions
# so it treats leaves outside of ball to see to which partition they belong to
def extender(graph, bipartitions, leaves,distance_matrix):
    u, v = leaves
    extended_bipartitions = []
    K = graph.copy()

    for bipartition in bipartitions:
        psi_u, psi_v = bipartition
        #print("psi_u, psi_v",  psi_u, psi_v)
        for node_u in psi_u:
            for node_v in psi_v:
                #print(node_u, node_v)
                if K.has_edge(node_u, node_v):
                    #print("remove edge between ", node_u, node_v)
                    K.remove_edge(node_u, node_v)

        #Initialize the extended partition 
        extended_psi_u = set(psi_u)
        extended_psi_v = set(psi_v)

        #For each node w not in psi_u union psi_v, add w to the partition it's connected to in K
        #print("set(K.nodes) - (psi_u union psi_v)", set(K.nodes) - (psi_u | psi_v))
        for w in set(K.nodes) - (psi_u | psi_v):
            max_dist_to_psi_u = max(distance_matrix[w][node_u] for node_u in psi_u)
            max_dist_to_psi_v = max(distance_matrix[w][node_v] for node_v in psi_v)
            #print("max_dist_to_psi: ", max_dist_to_psi_u, max_dist_to_psi_v)
            # Add w to the partition it's closer to, based on the distance m
            if max_dist_to_psi_v >= max_dist_to_psi_u:
                extended_psi_v.add(w)
            #print(w, " added to bipartition v ", extended_psi_v)
            else:
                extended_psi_u.add(w)
            #print(w, " added to bipartition u ", extended_psi_u)

        # Each side of the partition should now be extended
        extended_bipartitions.append((extended_psi_u, extended_psi_v))
    return extended_bipartitions


# random.shuffle(unique_sorted_tuples)
# print("reconstructed unrooted tree")
# print(reconstruct_unrooted_from_bipartitions(unique_sorted_tuples, Ntips).as_string(schema="newick"))
# recovered_tree_string = reconstruct_unrooted_from_bipartitions(unique_sorted_tuples, Ntips).as_string(schema="newick")
# recovered_newick_tree = re.search(r' (?=\()(.+;)', recovered_tree_string).group(1)
# recovered_newick_tree

#get unique trees from all tuples, shuffles bipartitions "runs" times to construct a tree (greedy manner, so need to shuffle)
def get_unique_trees(unique_sorted_tuples, Ntips, runs=10):
    unique_trees = set()
    taxon_namespace = dendropy.TaxonNamespace()

    for _ in range(runs):
        random.shuffle(unique_sorted_tuples)

        recovered_tree = reconstruct_unrooted_from_bipartitions(unique_sorted_tuples, Ntips)
        #print("tree: ", recovered_tree)
        #print('colless: ', treemeasure.colless_tree_imbalance(recovered_tree))
        recovered_tree_string = recovered_tree.as_string(schema="newick")
        recovered_newick_tree = re.search(r'(?=\()(.+;)', recovered_tree_string).group(1)

        is_unique = True
        for existing_tree_newick in unique_trees:
            existing_tree = dendropy.Tree.get(data=existing_tree_newick, schema="newick", taxon_namespace=taxon_namespace)
            new_tree = dendropy.Tree.get(data=recovered_newick_tree, schema="newick", taxon_namespace=taxon_namespace)
            rf_distance = dendropy.calculate.treecompare.symmetric_difference(existing_tree, new_tree)
            if rf_distance == 0:
                is_unique = False
                break

        if is_unique:
            unique_trees.add(recovered_newick_tree)

    return unique_trees


def get_internal_edge_lengths(tree):
    internal_edge_lengths = []
    for node in tree.traverse():
        if not node.is_leaf() and node.dist is not None:
        #if node.dist is not None:
            internal_edge_lengths.append(node.dist)
    return internal_edge_lengths




def jc69_distance(sequence1, sequence2):
    n = len(sequence1)
    diffs = sum(1 for x, y in zip(sequence1, sequence2) if x != y)

    p = diffs / n

    if p < 0.75:
        distance = -3/4 * np.log(1 - 4/3 * p)
    else:
        distance = float('inf')  

    return distance

def calculate_jc69_distance_matrix(msa):
    num_seqs = len(msa)
    labels = [record.id for record in msa]
    matrix = []
    for i in range(num_seqs):
        row = []
        for j in range(i + 1):  
            if i == j:
                row.append(0.0)  
            else:
                dist = jc69_distance(str(msa[j].seq), str(msa[i].seq))
                row.append(dist)
        matrix.append(row)
    return DistanceMatrix(names=labels, matrix=matrix)



def is_multiburcating_unrooted_tree(newick_string):
    tree = Phylo.read(StringIO(newick_string), "newick")
    is_multiburcating = any(len(node.clades) > 3 for node in tree.find_clades())

    return is_multiburcating


def extract_numbers(s):
    return re.findall(r'\d+', s)


def create_mappings(unique_sorted_tuples):
    unique_numbers = set()
    for A, B in unique_sorted_tuples:
        unique_numbers.update(A, B)
    sorted_numbers = sorted(unique_numbers)
    # Mapping from original to continuous range
    mapping = {number: i for i, number in enumerate(sorted_numbers)}
    
    # Reverse mapping from continuous range back to original
    reverse_mapping = {i: number for number, i in mapping.items()}

    return mapping, reverse_mapping

def find_identical_species(distance_matrix):
    identical_species = False
    for i in range(distance_matrix.shape[0]):
        for j in range(i + 1, distance_matrix.shape[1]):
            if distance_matrix.iloc[i, j] == 0:
                print(f"Species {distance_matrix.index[i]} and {distance_matrix.columns[j]} are identical.")
                identical_species = True
    if not identical_species:
        print("No identical species found.")
        return 1
    else:
        return 0
        
#group results by tau
def process_forest_results_grouped_by_tau(results):
    results_by_tau = {}
    for result in results:
        tau = result['tau']
        if tau not in results_by_tau:
            results_by_tau[tau] = []
        results_by_tau[tau].append(result)

    processed_results = {}
    for tau, group in results_by_tau.items():
        # Sort group by M in descending order
        sorted_group = sorted(group, key=lambda x: -x['M'])

        # find an entry with all induced_rf values are 0
        selected_result = None
        for result in sorted_group:
            if all(rf == 0 for rf in (result['induced_rfs'] if isinstance(result['induced_rfs'], list) else [result['induced_rfs']])):
                selected_result = result
                break

        # If not found, use the one with the largest M in the group
        if selected_result is None:
            selected_result = sorted_group[0]

        processed_results[tau] = selected_result

    return processed_results


def calculate_chord_depth(newick_string, distance_matrix):
    # Load the true tree from the Newick string
    tree = dendropy.Tree.get(data=newick_string, schema="newick")
    max_min_path_length = 0

    # function to collect leaf node labels
    def collect_leaf_labels(node):
        return [leaf.taxon.label for leaf in node.leaf_iter()]

    # Iterate over all edges in the tree
    for edge in tree.postorder_edge_iter():
        if not edge.is_terminal():
            min_path_length = float('inf')
            
            # Temporarily remove the edge and get the two components
            child_node = edge.head_node
            parent_node = edge.tail_node
            if parent_node and child_node:
                parent_node.remove_child(child_node)
                
                component1 = set(collect_leaf_labels(child_node))
                component2 = set(collect_leaf_labels(tree.seed_node)) - component1  # seed_node is the root node of the tree

                # Restore the edge
                parent_node.add_child(child_node)
                
                # Compute distances between all pairs from different components
                for taxon1 in component1:
                    for taxon2 in component2:
                        try:
                            # Safely access the distance from the matrix
                            path_length = distance_matrix[taxon1][taxon2]
                            # Update the minimum path length for this edge
                            if path_length < min_path_length:
                                min_path_length = path_length
                        except KeyError:
                            print(f"Missing distance data for taxa: {taxon1}, {taxon2}")

            # Update the overall maximum of the minimum path lengths
            if min_path_length < float('inf') and min_path_length > max_min_path_length:
                max_min_path_length = min_path_length
    
    return max_min_path_length

#below is the test to compute the chord depth
# newick_string = "((A,B),(C,D));"
# distance_matrix = {
#     'A': {'A': 0, 'B': 5, 'C': 19, 'D': 9},
#     'B': {'A': 5, 'B': 0, 'C': 10, 'D': 10},
#     'C': {'A': 19, 'B': 10, 'C': 0, 'D': 3},
#     'D': {'A': 9, 'B': 10, 'C': 3, 'D': 0}
# }
# chord_depth = calculate_chord_depth(newick_string, distance_matrix)
# print("Chord Depth of the tree:", chord_depth)
# calculate_chord_depth(true_tree_string, true_distance_matrix)


def prune_tree_and_remove_branch_lengths(tree, species_list):
    #tree = Phylo.read(StringIO(newick), 'newick')
    tree.prune(species_list, preserve_branch_length=False)
    return tree.write(format=9)

#refer to Daskalakis, Mossel, Roth paper (2011, Prune Deep...)
def get_forest_m(true_newick_string, nj_tree, true_a, a, m_steps, taus, M_steps, Ntips, Thr, k, epsilon, Nruns, taxon_namespace):
    collected_data = []
    t = Tree(true_newick_string)
    njt = Tree(nj_tree)
    if find_identical_species(a) == 1:
        M_min = a.min().min() + epsilon
        M_max = a.max().max() + epsilon
        Ms = np.linspace(M_min, M_max, M_steps) 
        ms = np.linspace(M_min, M_max, m_steps) 
        for m in ms:
            for tau in taus:
                for M in Ms:
                    if M > 2*m + 3*tau: #Theorem 1 condition
                        true_g = construct_clustering_graph(np.array(a), m)
                        num_components = nx.number_connected_components(true_g)
                        components = list(nx.connected_components(true_g))
                        #print(f"Processing with m = {m}, tau = {tau} and M = {M}")
                        if all(len(component) > Thr  for component in components): #Thr = 2
                            all_components_single_tree = True
                            component_trees = {}
                            component_sizes = []
                            for i, component in enumerate(components):
                                component_subgraph = true_g.subgraph(component)
                                component_size = len(component_subgraph.nodes())
                                component_sizes.append(component_size)
                                leaves_pairs = [(u, v) for u in component_subgraph.nodes() for v in component_subgraph.nodes() if u < v]
                                final_bipartitions = []
                                mini_bipartitions_list = []
                                
                                for leaves in leaves_pairs:
                                    mini_bipartitions, clusters = mini_contractor(component_subgraph, np.array(a),  leaves, M, tau)
                                    #print("mb: ", M, mini_bipartitions)
                                    mini_bipartitions_list.append(mini_bipartitions)
                                    extended_bipartition = extender(component_subgraph, mini_bipartitions, leaves, np.array(a))
                                    #print("eb: ", tau, extended_bipartition)
                                    final_bipartitions.append(extended_bipartition)
                                flattened_bipartitions = [tuple(frozenset(part) for part in bipartition) for sublist in final_bipartitions for bipartition in sublist]
                                unique_bipartitions = list(set(flattened_bipartitions))
                                unique_bipartitions_as_sets = [tuple(set(part) for part in bipartition) for bipartition in unique_bipartitions]
                                unique_bipartitions_as_sets.sort()
                                sorted_frozenset_data = [sort_and_freeze(pair) for pair in unique_bipartitions_as_sets]
                                unique_sorted_frozenset_data = set(sorted_frozenset_data)
                                unique_sorted_tuples = [tuple(map(set, pair)) for pair in unique_sorted_frozenset_data]
                                if len(unique_sorted_tuples) > 0:
                                    unique_trees = get_unique_trees(unique_sorted_tuples, len(component_subgraph.nodes()), 10)
                                    if len(unique_trees) == 1:
                                        component_trees[i] = unique_trees
                                    else:
                                        all_components_single_tree = False
                                        break
                                else:
                                    all_components_single_tree = False
                                    break

                            if all_components_single_tree:
                                induced_rfs = []
                                tree_parts = []
                                induced_rfs_nj = []
                                nj_parts = []
                                #print(f"Combination M = {M}, tau = {tau} results in single unique trees for all components")
                                #print(f'diff < tau ', np.abs(true_a - a).max().max(), tau)
                                #print(f'd < M + tau ', true_a.max().max(), M + tau)
                                #print(f'd-hat < M + tau ', a.max().max(), M + tau)
                                for component_idx, trees in component_trees.items():
                                    for tree in trees:
                                        print(tree)  
                                        forest_part_tree = dendropy.Tree.get(data=tree, schema="newick", taxon_namespace=taxon_namespace)
                                        tree_parts.append(tree)
                                        species_list = extract_numbers(tree)
#                                            dist_part_obs_a = a.loc[species_list, species_list]
#                                            dist_part_true_a = true_a.loc[species_list, species_list]
#                                             diff = np.abs(dist_part_true_a - dist_part_obs_a)
#                                             print(f'diff  tau ', diff.max().max(), tau)
#                                             if diff.max().max() < tau:
#                                                 print(f'diff < tau ', diff.max().max(), tau)
#                                                 print(f'd < M + tau ', dist_part_true_a.max().max(), M + tau)
#                                                 print(f'd-hat < M + tau ', dist_part_obs_a.max().max(), M + tau)
#                                             else: 
#                                                 print('not a distortion: ', diff.max().max())
                                        true_tree = t.copy()
                                        #print("species_list ", species_list)
                                        true_tree.prune(species_list, preserve_branch_length=False)
                                        #print(true_tree.write())
                                        true_part_tree = dendropy.Tree.get(data=true_tree.write(), schema="newick", taxon_namespace=taxon_namespace)
                                        induced_rf = dendropy.calculate.treecompare.symmetric_difference(true_part_tree, forest_part_tree)
                                        induced_rfs.append(induced_rf)
                                        nj_t = njt.copy()
                                        nj_t = prune_tree_and_remove_branch_lengths(nj_t, species_list)
                                        nj_parts.append(nj_t)
                                        nj_part_tree = dendropy.Tree.get(data=nj_t, schema="newick", taxon_namespace=taxon_namespace)
                                        induced_rf_nj = dendropy.calculate.treecompare.symmetric_difference(true_part_tree, nj_part_tree)
                                        induced_rfs_nj.append(induced_rf_nj)    
                                #print(f"pd_rf: {induced_rf}, nj_rf: {induced_rf_nj}")
                                component_info = {
                                    'M': M,
                                    'm': m,
                                    'tau': tau,
                                    'max_dhat': a.max().max(),
                                    'max_d': true_a.max().max(),
                                    'either_d_dhat_less_Mplustau': 1 if true_a.max().max() < M + tau or a.max().max() < M + tau else 0,
                                    'd_minus_dhat_less_tau': 1 if np.abs(true_a - a).max().max() < tau else 0,
                                    'num_components': num_components,
                                    'component_sizes': component_sizes if num_components > 1 else component_sizes[0],
                                    'induced_rfs': induced_rfs if num_components > 1 else induced_rfs[0],
                                    'forest': tree_parts if num_components > 1 else tree_parts[0],
                                    'induced_rfs_nj': induced_rfs_nj if num_components > 1 else induced_rfs_nj[0],
                                    'nj_parts': nj_parts if num_components > 1 else nj_parts[0]
                                }
                                if num_components > 1:
                                    #print(f"Number of components: {num_components}, m: {m}, M: {M}")
                                    #print(component_info)
                                    collected_data.append(component_info)

                #             else:
                #                 print(f"Combination M = {M}, tau = {tau} does not satisfy the condition")
                
        else:
            print("Identical species found.")
    return collected_data

def get_NJ(true_newick_string, a, Ntips, taxon_namespace):
    number_strings = [str(i) for i in range(Ntips)]
    full_matrix = a.values.tolist()
    #taxa = a.index.tolist()  
    lower_triangular_matrix = [row[:i+1] for i, row in enumerate(full_matrix)]
    dm = DistanceMatrix(names=number_strings, matrix=lower_triangular_matrix)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    nj_tree = io.StringIO()
    Phylo.write(tree, nj_tree, "newick")
    nj_tree = re.sub(r'Inner\d+:', '', nj_tree.getvalue())
    #print(f'NJ tree: {nj_tree}')
    nj_newick_tree = dendropy.Tree.get(data=nj_tree, schema="newick", taxon_namespace=taxon_namespace)
    true_newick_tree = dendropy.Tree.get(data=true_newick_string, schema="newick", taxon_namespace=taxon_namespace)
    return nj_tree, dendropy.calculate.treecompare.symmetric_difference(true_newick_tree, nj_newick_tree)




def get_forest_pars(true_newick_string, nj_tree, true_a, a, m_steps, taus, M_steps, Ntips, Thr, k, epsilon, Nruns, taxon_namespace):
    collected_data = []
    t = Tree(true_newick_string)
    njt = Tree(nj_tree)
    if find_identical_species(a) == 1:
        chord_d = calculate_chord_depth(true_newick_string, true_a)
        M_min = chord_d/2 - 0.1
        M_max  = chord_d + 1
        #M_min = a.min().min() + epsilon
        #M_max = a.max().max() + epsilon
        Ms = np.linspace(M_min, M_max, M_steps) 
        ms = np.linspace(M_min, M_max, m_steps) 
        print(f'chord_true_depth: {chord_d}, Ms: {M_min}, {M_max}')
        
        #chord_d = calculate_chord_depth(true_newick_string, true_a)
#         M_min = a.min().min() + epsilon
#         M_max = a.max().max() + epsilon
#         Ms = np.linspace(M_min, M_max, M_steps) 
#         ms = np.linspace(M_min, M_max, m_steps) 
#         print(f'chord_est. Ms: {M_min}, {M_max}')
    
        for m in ms:
            for tau in taus:
                for M in Ms:
                    if M > 2*m + 3*tau:
                        true_g = construct_clustering_graph(np.array(a), m)
                        num_components = nx.number_connected_components(true_g)
                        #print(f'num_components: {num_components}, m: {m}')
                        components = list(nx.connected_components(true_g))
                        #print(f"Processing with m = {m}, tau = {tau} and M = {M} and ChordD: {chord_d}")
                        if all(len(component) > Thr  for component in components): #Thr = 2
                            #print(f"Processing with m = {m}, tau = {tau} and M = {M} and ChordD: {chord_d}")
                            all_components_single_tree = True
                            component_trees = {}
                            component_sizes = []
                            for i, component in enumerate(components):
                                component_subgraph = true_g.subgraph(component)
                                component_size = len(component_subgraph.nodes())
                                component_sizes.append(component_size)
#                                 print(f"component_size {component_size}")
#                                 plt.figure(figsize=(4, 4))
#                                 nx.draw(component_subgraph, with_labels=True, font_weight='bold')
#                                 plt.show()
                                
                                leaves_pairs = [(u, v) for u in component_subgraph.nodes() for v in component_subgraph.nodes() if u < v]

                                final_bipartitions = []
                                mini_bipartitions_list = []

                                for leaves in leaves_pairs:
                                    mini_bipartitions, clusters = mini_contractor(component_subgraph, np.array(a),  leaves, M, tau)
                                    #print("mb: ", mini_bipartitions)
                                    mini_bipartitions_list.append(mini_bipartitions)
                                    extended_bipartition = extender(component_subgraph, mini_bipartitions, leaves, np.array(a))
                                    #print("eb: ", extended_bipartition)
                                    final_bipartitions.append(extended_bipartition)
                                flattened_bipartitions = [tuple(frozenset(part) for part in bipartition) for sublist in final_bipartitions for bipartition in sublist]
                                unique_bipartitions = list(set(flattened_bipartitions))
                                unique_bipartitions_as_sets = [tuple(set(part) for part in bipartition) for bipartition in unique_bipartitions]
                                unique_bipartitions_as_sets.sort()
                                sorted_frozenset_data = [sort_and_freeze(pair) for pair in unique_bipartitions_as_sets]
                                unique_sorted_frozenset_data = set(sorted_frozenset_data)
                                unique_sorted_tuples = [tuple(map(set, pair)) for pair in unique_sorted_frozenset_data]
                                if len(unique_sorted_tuples) > 0:
                                    unique_trees = get_unique_trees(unique_sorted_tuples, len(component_subgraph.nodes()), 10)
                                    if len(unique_trees) == 1:
                                        component_trees[i] = unique_trees
                                    else:
                                        all_components_single_tree = False
                                        #print(f"Processing with m = {m}, tau = {tau} and M = {M} and ChordD: {chord_d}")
                                        #print("Multiple trees from Tree Popping")
                                        #print(unique_trees)
                                        break
                                else:
                                    all_components_single_tree = False
                                    break

                            if all_components_single_tree:
                                induced_rfs = []
                                tree_parts = []
                                induced_rfs_nj = []
                                nj_parts = []
                                #print(f"Combination M = {M}, m = {m}, tau = {tau} results in single unique trees for all components")

                                #print(f'diff < tau ', np.abs(true_a - a).max().max(), tau)
                                #print(f'd < M + tau ', true_a.max().max(), M + tau)
                                #print(f'd-hat < M + tau ', a.max().max(), M + tau)

                                for component_idx, trees in component_trees.items():
#                                     print(f"Num comps {num_components}")
#                                     print(f"Component {component_idx + 1}, size {component_size}")
#                                     plt.figure(figsize=(8, 8))
#                                     nx.draw(true_g, with_labels=True, font_weight='bold')
#                                     plt.show()
                                    for tree in trees:
                                        #print(tree)  
                                        forest_part_tree = dendropy.Tree.get(data=tree, schema="newick", taxon_namespace=taxon_namespace)
                                        tree_parts.append(tree)
                                        species_list = extract_numbers(tree)
#                                        dist_part_obs_a = a.loc[species_list, species_list]
#                                        dist_part_true_a = true_a.loc[species_list, species_list]
#                                             diff = np.abs(dist_part_true_a - dist_part_obs_a)
#                                             print(f'diff  tau ', diff.max().max(), tau)
#                                             if diff.max().max() < tau:
#                                                 print(f'diff < tau ', diff.max().max(), tau)
#                                                 print(f'd < M + tau ', dist_part_true_a.max().max(), M + tau)
#                                                 print(f'd-hat < M + tau ', dist_part_obs_a.max().max(), M + tau)
#                                             else: 
#                                                 print('not a distortion: ', diff.max().max())
                                        true_tree = t.copy()
                                        #print("species_list ", species_list)
                                        true_tree.prune(species_list, preserve_branch_length=False)
                                        #print(true_tree.write())
                                        true_part_tree = dendropy.Tree.get(data=true_tree.write(), schema="newick", taxon_namespace=taxon_namespace)
                                        induced_rf = dendropy.calculate.treecompare.symmetric_difference(true_part_tree, forest_part_tree)
                                        induced_rfs.append(induced_rf)
                                        nj_t = njt.copy()
                                        nj_t = prune_tree_and_remove_branch_lengths(nj_t, species_list)
                                        nj_parts.append(nj_t)

                                        nj_part_tree = dendropy.Tree.get(data=nj_t, schema="newick", taxon_namespace=taxon_namespace)
                                        induced_rf_nj = dendropy.calculate.treecompare.symmetric_difference(true_part_tree, nj_part_tree)
                                        induced_rfs_nj.append(induced_rf_nj)    
                                #print(f"pd_rf: {induced_rf}, nj_rf: {induced_rf_nj}")
                                component_info = {
                                    'M': M,
                                    'm': m,
                                    'MmRange': [M_min, M_max],
                                    'tau': tau,
                                    'max_dhat': a.max().max(),
                                    'max_d': true_a.max().max(),
                                    'chord_true_depth': chord_d,
                                    'either_d_dhat_less_Mplustau': 1 if true_a.max().max() < M + tau or a.max().max() < M + tau else 0,
                                    'd_minus_dhat_less_tau': 1 if np.abs(true_a - a).max().max() < tau else 0,
                                    'num_components': num_components,
                                    'component_sizes': component_sizes if num_components > 1 else component_sizes[0],
                                    'induced_rfs': induced_rfs if num_components > 1 else [induced_rfs[0]],
                                    'forest': tree_parts if num_components > 1 else tree_parts[0],
                                    'induced_rfs_nj': induced_rfs_nj if num_components > 1 else [induced_rfs_nj[0]],
                                    'nj_parts': nj_parts if num_components > 1 else nj_parts[0]
                                }
                                print(component_info)
                                #if num_components > 1:
                                    #print(f"Number of components: {num_components}, m: {m}, M: {M}")
                                    #print(component_info)
                                collected_data.append(component_info)
                            else:
                                component_info = {
                                    'M': M,
                                    'm': m,
                                    'MmRange': [M_min, M_max],
                                    'tau': tau,
                                    'max_dhat': a.max().max(),
                                    'max_d': true_a.max().max(),
                                    'chord_true_depth': chord_d,
                                    'either_d_dhat_less_Mplustau': 1 if true_a.max().max() < M + tau or a.max().max() < M + tau else 0,
                                    'd_minus_dhat_less_tau': 1 if np.abs(true_a - a).max().max() < tau else 0,
                                    'num_components': num_components,
                                    'component_sizes': component_sizes if num_components > 1 else component_sizes[0],
                                    'induced_rfs': None,
                                    'forest': None,
                                    'induced_rfs_nj': None,
                                    'nj_parts': None
                                }
                                collected_data.append(component_info)
                #                 print(f"Combination M = {M}, tau = {tau} does not satisfy the condition")
                            

    else:
        print("Identical species found.")

    return collected_data



def adjust_taxon_labels(newick):
    def replace_label(match):
        num = int(match.group(1))  
        return str(num - 1)  
    adjusted_newick = re.sub(r'T(\d+)', replace_label, newick)
    return adjusted_newick


def adjust_labels(labels):
    new_labels = []
    for label in labels:
        if label.startswith('T') and label[1:].isdigit():
            new_label = str(int(label[1:]) - 1)
        else:
            new_label = label  # Keep the label unchanged if it doesn't match expected format
        new_labels.append(new_label)
    return new_labels


def get_forest_custom_input(true_newick_string, nj_tree, true_a, a, ms, taus, Ms, Ntips, Thr, k, epsilon, Nruns, taxon_namespace):
    collected_data = []
    t = Tree(true_newick_string)
    njt = Tree(nj_tree)
    chord_d = calculate_chord_depth(true_newick_string, true_a)
    #print(f'chord depth: {chord_d}')
    if find_identical_species(a) == 1:    
        for m in ms:
            for tau in taus:
                for M in Ms:
                    if M > 2*m + 3*tau:
                        true_g = construct_clustering_graph(np.array(a), m)
                        num_components = nx.number_connected_components(true_g)
                        components = list(nx.connected_components(true_g))
                        if all(len(component) > Thr  for component in components): #Thr = 2
                            all_components_single_tree = True
                            component_trees = {}
                            component_sizes = []
                            for i, component in enumerate(components):
                                component_subgraph = true_g.subgraph(component)
                                component_size = len(component_subgraph.nodes())
                                component_sizes.append(component_size)
                                leaves_pairs = [(u, v) for u in component_subgraph.nodes() for v in component_subgraph.nodes() if u < v]

                                final_bipartitions = []
                                mini_bipartitions_list = []

                                for leaves in leaves_pairs:
                                    #print(f'leaves: {leaves}')
                                    mini_bipartitions, clusters = mini_contractor(component_subgraph, np.array(a),  leaves, M, tau)
                                    #print("mb: ", mini_bipartitions)
                                    mini_bipartitions_list.append(mini_bipartitions)
                                    extended_bipartition = extender(component_subgraph, mini_bipartitions, leaves, np.array(a))
                                    #print("eb: ", extended_bipartition)
                                    final_bipartitions.append(extended_bipartition)
                                flattened_bipartitions = [tuple(frozenset(part) for part in bipartition) for sublist in final_bipartitions for bipartition in sublist]
                                unique_bipartitions = list(set(flattened_bipartitions))
                                unique_bipartitions_as_sets = [tuple(set(part) for part in bipartition) for bipartition in unique_bipartitions]
                                unique_bipartitions_as_sets.sort()
                                sorted_frozenset_data = [sort_and_freeze(pair) for pair in unique_bipartitions_as_sets]
                                unique_sorted_frozenset_data = set(sorted_frozenset_data)
                                unique_sorted_tuples = [tuple(map(set, pair)) for pair in unique_sorted_frozenset_data]
                                if len(unique_sorted_tuples) > 0:
                                    unique_trees = get_unique_trees(unique_sorted_tuples, len(component_subgraph.nodes()), 10)
                                    if len(unique_trees) == 1:
                                        component_trees[i] = unique_trees
                                    else:
                                        all_components_single_tree = False
                                        break
                                else:
                                    all_components_single_tree = False
                                    break

                            if all_components_single_tree:
                                induced_rfs = []
                                tree_parts = []
                                induced_rfs_nj = []
                                nj_parts = []

                                for component_idx, trees in component_trees.items():
                                    for tree in trees:
                                        forest_part_tree = dendropy.Tree.get(data=tree, schema="newick", taxon_namespace=taxon_namespace)
                                        tree_parts.append(tree)
                                        species_list = extract_numbers(tree)
                                       # dist_part_obs_a = a.loc[species_list, species_list]
                                       # dist_part_true_a = true_a.loc[species_list, species_list]
                                        true_tree = t.copy()
                                        true_tree.prune(species_list, preserve_branch_length=False)
                                        true_part_tree = dendropy.Tree.get(data=true_tree.write(), schema="newick", taxon_namespace=taxon_namespace)
                                        induced_rf = dendropy.calculate.treecompare.symmetric_difference(true_part_tree, forest_part_tree)
                                        induced_rfs.append(induced_rf)
                                        nj_t = njt.copy()
                                        nj_t = prune_tree_and_remove_branch_lengths(nj_t, species_list)
                                        nj_parts.append(nj_t)

                                        nj_part_tree = dendropy.Tree.get(data=nj_t, schema="newick", taxon_namespace=taxon_namespace)
                                        induced_rf_nj = dendropy.calculate.treecompare.symmetric_difference(true_part_tree, nj_part_tree)
                                        induced_rfs_nj.append(induced_rf_nj)    
                                component_info = {
                                    'M': M,
                                    'm': m,
                                    'tau': tau,
                                    'max_dhat': a.max().max(),
                                    'max_d': true_a.max().max(),
                                    'chord_true_depth': chord_d,
                                    'd_less_Mplustau': 1 if true_a.max().max() < M + tau else 0,
                                    'd_dhat_less_Mplustau': 1 if a.max().max() < M + tau else 0,
                                    'd_minus_dhat_less_tau': 1 if np.abs(true_a - a).max().max() < tau else 0,
                                    'num_components': num_components,
                                    'component_sizes': component_sizes if num_components > 1 else component_sizes[0],
                                    'induced_rfs': induced_rfs if num_components > 1 else [induced_rfs[0]],
                                    'forest': tree_parts if num_components > 1 else tree_parts[0],
                                    'induced_rfs_nj': induced_rfs_nj if num_components > 1 else [induced_rfs_nj[0]],
                                    'nj_parts': nj_parts if num_components > 1 else nj_parts[0]
                                }
                                print(component_info)
                                collected_data.append(component_info)
                            else:
                                component_info = {
                                    'M': M,
                                    'm': m,
                                    'tau': tau,
                                    'max_dhat': a.max().max(),
                                    'max_d': true_a.max().max(),
                                    'chord_true_depth': chord_d,
                                    'd_less_Mplustau': 1 if true_a.max().max() < M + tau else 0,
                                    'd_dhat_less_Mplustau': 1 if a.max().max() < M + tau else 0,                                    'd_minus_dhat_less_tau': 1 if np.abs(true_a - a).max().max() < tau else 0,
                                    'num_components': num_components,
                                    'component_sizes': component_sizes if num_components > 1 else component_sizes[0],
                                    'induced_rfs': None,
                                    'forest': None,
                                    'induced_rfs_nj': None,
                                    'nj_parts': None
                                }
                                collected_data.append(component_info)                            
    else:
        print("Identical species found.")

    return collected_data
