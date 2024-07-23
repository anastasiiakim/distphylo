#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run forest output, needs IQTREE2

"""
from forest_functions import *

#True dist example  m < chord depth, PD outperforms NJ
np.random.seed(2024)


Ntips = [32]
lb = 0.05
ub = 0.1
Thr = 2
taus = np.linspace(0.02, lb/2, 4)
epsilon = 1e-8
Nruns = 5
k_values = [64]
num_aln_reps = 2 #should be greater than 1 or otherwise IQTREE will not produce correct file name
forest_results = {}
Ms = [1.5, 1, 0.8]
ms = [0.6, 0.5, 0.4, 0.3, 0.25]
ntips = Ntips[0]

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

tree_path = os.path.join(tipdir, f"true_tree_{ntips}.tre")
with open(tree_path, 'w') as tree_file:
    tree_file.write(true_tree.as_string(schema="newick"))

draw_true_tree = Phylo.read(tree_path, 'newick')
Phylo.draw(draw_true_tree, branch_labels = lambda c: np.round(c.branch_length, 3))

taxon_namespace = dendropy.TaxonNamespace()
true_tree_string = re.search(r' (?=\()(.+;)', true_tree.as_string(schema="newick")).group(1) + '\n'
true_tree_string = adjust_taxon_labels(true_tree_string)
internal_edges = get_internal_edge_lengths(Tree(true_tree_string))
print(f"internal_edges: {np.sort(internal_edges)}")
print(f"0.5*internal_edges: {0.5*np.sort(internal_edges)}")
print(f'chord_true_depth: {calculate_chord_depth(true_newick_string, true_distance_matrix)}')
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
#print(true_distance_matrix)
true_distance_matrix.index = true_distance_matrix.index.astype(int)  # Convert index to integer
true_distance_matrix = true_distance_matrix.sort_index()  # Sort DataFrame by index
true_distance_matrix.index = true_distance_matrix.index.astype(str)  

true_distance_matrix.columns = true_distance_matrix.columns.astype(int)  # Convert column names to integer
true_distance_matrix = true_distance_matrix.sort_index(axis=1)  # Sort DataFrame by column names
true_distance_matrix.columns = true_distance_matrix.columns.astype(str)  # Convert column names back to string if necessary
#print(true_distance_matrix)


iqtree2_path = "/Users/akim/Documents/phylomc/iqtree2"

treefile = os.path.join(tipdir, f"true_tree_{ntips}.tre")

for k in k_values:
    print(f"seq_length: {k}")
    alnfile = os.path.join(tipdir, f"aln_{k}")
    cmd = [iqtree2_path, "--alisim", alnfile, "-t", treefile, "-m", "JC", "--length", str(k), "--num-alignments", str(num_aln_reps), "--out-format", "fasta", "--seed", "0"]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    for i_replicate in range(2, num_aln_reps+1):
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

        print(f'diff < tau ', np.abs(true_distance_matrix - est_dist_mat).max().max())
        print(f'd < M + tau ', true_distance_matrix.max().max())
        print(f'd-hat < M + tau ', est_dist_mat.max().max())


        constructor = DistanceTreeConstructor()
        tree = constructor.nj(est_dist_mat_orig)
        nj_tree = io.StringIO()
        Phylo.write(tree, nj_tree, "newick") 
        nj_tree_string = re.sub(r'Inner\d+:', '', nj_tree.getvalue())
        nj_tree_string = adjust_taxon_labels(nj_tree_string)
        nj_tree = dendropy.Tree.get(data=nj_tree_string.rstrip('\n'), schema="newick", taxon_namespace=taxon_namespace)
        rf_distance = dendropy.calculate.treecompare.symmetric_difference(true_tree, nj_tree)
        print(f'rf_distance: {rf_distance}')
#             key = f"forest_results_{i_replicate}"
#             forest_results[key] = get_forest_pars(true_tree_string, nj_tree_string, true_distance_matrix, est_dist_mat, m_steps, taus, M_steps, ntips, Thr, k, epsilon, Nruns)
        key = f"forest_results_{i_replicate}"
        Ms = [1.5, 1, 0.8]
        ms = [0.6, 0.5, 0.4, 0.3, 0.25]
        forest_results[key] = get_forest_custom_input(true_tree_string, nj_tree_string, true_distance_matrix, true_distance_matrix, ms, taus, Ms, ntips, Thr, k, epsilon, Nruns, printout = 0)


        #Est dist example  m < chord depth
os.chdir("/Users/akim/Documents/phylomc/")
np.random.seed(2024)


Ntips = [32]
lb = 0.05
ub = 0.1
Thr = 2
taus = np.linspace(0.01, lb/2 + 0.03, 10)
epsilon = 1e-8
Nruns = 5
k_values = [512]
num_aln_reps = 2 #should be greater than 1 or otherwise IQTREE will not produce correct file name
forest_results = {}

ntips = Ntips[0]

tipdir = f"./rand_nj_lowestrate/{ntips}tips"
shutil.rmtree(tipdir)
os.makedirs(tipdir)
print(tipdir)
true_tree = dendropy.simulate.treesim.birth_death_tree(birth_rate=1.0, death_rate=0, num_extant_tips=ntips)
for edge in true_tree.edges():
   edge.length = np.random.uniform(lb, ub)
#print(true_tree)
for idx, leaf in enumerate(true_tree.leaf_nodes()):
   leaf.taxon.label = str(idx)

tree_path = os.path.join(tipdir, f"true_tree_{ntips}.tre")
with open(tree_path, 'w') as tree_file:
   tree_file.write(true_tree.as_string(schema="newick"))
draw_true_tree = Phylo.read(tree_path, 'newick')
Phylo.draw(draw_true_tree, branch_labels = lambda c: np.round(c.branch_length, 3))

taxon_namespace = dendropy.TaxonNamespace()
true_tree_string = re.search(r' (?=\()(.+;)', true_tree.as_string(schema="newick")).group(1) + '\n'
true_tree_string = adjust_taxon_labels(true_tree_string)



# true_tree_string = "(((T24:0.05219042818734324,(((T10:0.08397002617625708,T4:0.07369228517041093):0.08636200718422274,(T29:0.050955347393623306,((T30:0.09808887876540949,T7:0.08321843236782378):0.08012242695040633,T15:0.08033148096593382):0.08762991686042897):0.07241479122401523):0.055303143723813286,((T32:0.08350871484463479,T27:0.08678832962259567):0.06126770815963478,(T18:0.05477710769301827,T12:0.09804548711183389):0.06289978190390969):0.07245756574658582):0.060250947619714446):0.0594075980019253,(((T20:0.08989616985574918,T1:0.07720185992002081):0.08841269673313926,T16:0.0691353815323321):0.06410825597218217,((T28:0.08701340765703286,(T11:0.0718860852324932,T13:0.0941769351388293):0.061949341621786454):0.06429136942385916,(((T21:0.07088926924672168,T6:0.061288438379804154):0.08794768283867953,T8:0.07100490694652993):0.08922534285437106,(T9:0.07982163434604983,(T25:0.09462431931645277,T5:0.06002637219416135):0.09186618611554152):0.053218184544063574):0.06446405701663589):0.06908254751009987):0.06258836433840057):0.08495543738407912,((T19:0.06279604656941354,T14:0.09336161714831452):0.09476909222306378,((T3:0.07639526932429701,((T2:0.05320091880793926,(T22:0.09370199001871647,(T31:0.09998706533431981,T17:0.06734019851733619):0.05818336452764336):0.0951052352293099):0.06229742217631054,T23:0.06564390796720793):0.09616751934350093):0.07762484771855298,T26:0.09235520104952216):0.05082439674083155):0.0751197617182645):0.0794007259447699;"
# true_tree_string = adjust_taxon_labels(true_tree_string)
# taxon_namespace = dendropy.TaxonNamespace()
# true_tree = Tree.get(data=true_tree_string, schema="newick", taxon_namespace=taxon_namespace)
# draw_true_tree = Phylo.read(StringIO(true_tree_string), 'newick')
# Phylo.draw(draw_true_tree, branch_labels=lambda c: np.round(c.branch_length, 3))


internal_edges = get_internal_edge_lengths(Tree(true_tree_string))
print(f"internal_edges: {np.sort(internal_edges)}")
print(f"0.5*internal_edges: {0.5*np.sort(internal_edges)}")
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
#print(true_distance_matrix)
true_distance_matrix.index = true_distance_matrix.index.astype(int)  # Convert index to integer
true_distance_matrix = true_distance_matrix.sort_index()  # Sort DataFrame by index
true_distance_matrix.index = true_distance_matrix.index.astype(str)  

true_distance_matrix.columns = true_distance_matrix.columns.astype(int)  # Convert column names to integer
true_distance_matrix = true_distance_matrix.sort_index(axis=1)  # Sort DataFrame by column names
true_distance_matrix.columns = true_distance_matrix.columns.astype(str)  # Convert column names back to string if necessary
#print(true_distance_matrix)
print(f'chord_true_depth: {calculate_chord_depth(true_tree_string, true_distance_matrix)}')




iqtree2_path = "./qtree2"

treefile = os.path.join(tipdir, f"true_tree_{ntips}.tre")

for k in k_values:
    print(f"seq_length: {k}")
    alnfile = os.path.join(tipdir, f"aln_{k}")
    cmd = [iqtree2_path, "--alisim", alnfile, "-t", treefile, "-m", "JC", "--length", str(k), "--num-alignments", str(num_aln_reps), "--out-format", "fasta", "--seed", "0"]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    for i_replicate in range(2, num_aln_reps+1):
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

        print(f'diff < tau ', np.abs(true_distance_matrix - est_dist_mat).max().max())
        print(f'd < M + tau ', true_distance_matrix.max().max())
        print(f'd-hat < M + tau ', est_dist_mat.max().max())


        constructor = DistanceTreeConstructor()
        tree = constructor.nj(est_dist_mat_orig)
        nj_tree = io.StringIO()
        Phylo.write(tree, nj_tree, "newick") 
        nj_tree_string = re.sub(r'Inner\d+:', '', nj_tree.getvalue())
        nj_tree_string = adjust_taxon_labels(nj_tree_string)
        nj_tree = dendropy.Tree.get(data=nj_tree_string.rstrip('\n'), schema="newick", taxon_namespace=taxon_namespace)
        rf_distance = dendropy.calculate.treecompare.symmetric_difference(true_tree, nj_tree)
        print(f'rf_distance: {rf_distance}')
#             key = f"forest_results_{i_replicate}"
#             forest_results[key] = get_forest_pars(true_tree_string, nj_tree_string, true_distance_matrix, est_dist_mat, m_steps, taus, M_steps, ntips, Thr, k, epsilon, Nruns)
        key = f"forest_results_{i_replicate}"
        Ms = [1.5, 1.2, 1, 0.8]
        ms = [0.5, 0.4, 0.45, 0.3]#0.8, 0.7, 0.6, 0.5, 0.4
        forest_results[key] = get_forest_custom_input(true_tree_string, nj_tree_string, true_distance_matrix, est_dist_mat, ms, taus, Ms, ntips, Thr, k, epsilon, Nruns, printout = 0)
#line below is where mini contractor changed to not include v in the ball      
#get_forest_custom_input_changed(true_tree_string, nj_tree_string, true_distance_matrix, est_dist_mat, ms, taus, Ms, ntips, Thr, k, epsilon, Nruns, printout = 0)
#'(0,(14,(13,(11,12)),(20,(21,22),(19,(17,18))),(15,16)),(3,((9,10),(7,8)),(4,5,6),(1,2)));', '(23,24,(31,(25,(26,30,(27,28,29)))));'],


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

replicat = 2
negative_entries = [
    entry for entry in forest_results[f"forest_results_{replicat}"]
    if entry['induced_rfs'] is not None and any(rf is not None and rf != 0 for rf in entry['induced_rfs'])
]
positive_entries = [
    entry for entry in forest_results[f"forest_results_{replicat}"]
    if entry['induced_rfs'] is not None and all(rf is None or rf == 0 for rf in entry['induced_rfs'])
]

nonexistent_entries = [
    entry for entry in forest_results[f"forest_results_{replicat}"]
    if entry['induced_rfs'] is None
]
# Set up for 3D plot
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111, projection='3d')

# Plot positive entries if any
if positive_entries:
    xs = [entry['m'] for entry in positive_entries]
    ys = [entry['M'] for entry in positive_entries]
    zs = [entry['tau'] for entry in positive_entries]
    ax.scatter(xs, ys, zs, c='g', marker='o', label='No forest errors')

# Plot negative entries if any
if negative_entries:
    xs = [entry['m'] for entry in negative_entries]
    ys = [entry['M'] for entry in negative_entries]
    zs = [entry['tau'] for entry in negative_entries]
    ax.scatter(xs, ys, zs, c='b', marker='^', label='Forest errors or ambiguity')

# if nonexistent_entries:
#     xs = [entry['m'] for entry in nonexistent_entries]
#     ys = [entry['M'] for entry in nonexistent_entries]
#     zs = [entry['tau'] for entry in nonexistent_entries]
#     ax.scatter(xs, ys, zs, c='r', marker='^', label='NA')

# Adding labels and legend
#ax.set_xlabel('m: determines \n the forest size')
# ax.set_ylabel('M: determines long edges \n on the short paths')
# ax.set_zlabel(r'$\tau$: \n edges < $\tau$ are contracted')

ax.set_ylabel('')
ax.set_zlabel('')
ax.annotate('M: determines long edges\non the short paths',
            xy=(0.85, 0.1), xycoords='axes fraction',
            xytext=(10, 0), textcoords='offset points',
            ha='left', va='center', rotation=0)
ax.text2D(1.1, 0.5, r'$\tau$: edges < $\tau$'+'\n are contracted', transform=ax.transAxes, ha='left', rotation=0)
ax.text2D(0.15, -0.05, 'm: determines \n the forest size', transform=ax.transAxes, ha='left', rotation=0)

ax.legend(loc='upper right', bbox_to_anchor=(1.5, 1))
#plt.tight_layout(rect=[0, 0, 0.85, 1])


# Display the plot
plt.savefig('./32tips_mMtau_est_chord051.png', bbox_inches='tight', pad_inches=0.1)
plt.show()

#positive_entries
#forest_results[f"forest_results_{2}"]