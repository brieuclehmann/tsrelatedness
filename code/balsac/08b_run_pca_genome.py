import os
import sys
import time
import tskit
import msprime

import pandas as pd
import numpy as np

from scipy.sparse.linalg import LinearOperator,eigsh

sys.path.insert(0, "utils")
import pedigree_tools
import grm_tools

def genome_matvec(ts_list, A, x, mode):
    sample_x = A.T @ x
    return A @ sum(map(lambda ts: ts.genetic_relatedness_vector(sample_x, mode = mode, span_normalise = False), ts_list))

def ts_linop_pca(ts_list, mode='branch', n_pc=10):
    n_samples = ts_list[0].num_samples
    n_ind = n_samples // 2
    A = np.repeat(np.eye(n_ind), 2, axis = 1)
    #grm_matvec = lambda x: A @ sum(map(lambda ts: ts.genetic_relatedness_vector(x, mode = "branch", span_normalise = False), ts_list))
    grm_matvec = lambda x: genome_matvec(ts_list, A, x, mode)
    grm = LinearOperator((n_ind, n_ind), matvec=grm_matvec)
    eigval, eigvec = eigsh(grm, k=n_pc)
    prin_comp = (grm @ eigvec[:, ::-1])
    return prin_comp

def main():

    ###################
    #### Load data ####
    ###################

    pedigree_df = pd.read_csv("data/ascendance.csv", sep=";")
    region_df = pedigree_df[['FiveRegions','ProbandID']].drop_duplicates()
    region_df = region_df.rename(columns={'FiveRegions':'proband_region','ProbandID':'ind'})
    region_df = region_df.astype({"ind":'int'})
    full_df = pd.read_csv("data/balsac_pedigree.csv")
    proband_df = full_df[~(full_df['ind'].isin(full_df['mother'])) & ~(full_df['ind'].isin(full_df['father']))]
    proband_df = proband_df.merge(region_df, on='ind')

    # Extract metadata and sample nodes
    founder_df = pd.read_csv("data/balsac_founder_depths.csv")
    meta_df = region_df.merge(founder_df, left_on = 'ind', right_on = 'proband').drop(columns = 'proband')
    
    # Filter out probands with shallow average founder depth
    shallow_ind = founder_df[founder_df['average_depth'] < 3]['proband']
    proband_pca = list(proband_df['ind'])
    proband_pca = [i for i in proband_pca if i not in list(shallow_ind)]

    ####################
    #### branch PCA ####
    ####################
    ts_list = []
    override = True
    for chrom in range(1,23):
        ts_dir_out = f"output/balsac/chr{chrom}/sim1"
        ts = tskit.load(ts_dir_out + "/trees_recap.ts")
        ### Rescale tree sequence ###
        ts_path = ts_dir_out + "/trees_recap_rescaled.ts"
        if os.path.exists(ts_path) and not override:
            new_ts = tskit.load(ts_path)
        else:
            sample_nodes = [i.nodes for i in ts.individuals() if int(i.metadata['individual_name']) in proband_pca]
            sample_nodes = [node for pair in sample_nodes for node in pair]
            ts_sub = ts.simplify(sample_nodes)
            map_file_name = f'data/maps/genetic_map_GRCh37_chr{chrom}.txt'
            rate_map = msprime.intervals.RateMap.read_hapmap(map_file_name)
            t = ts_sub.dump_tables()
            t.mutations.clear()
            t.sites.clear()
            new_left = rate_map.get_cumulative_mass(t.edges.left)
            new_right = rate_map.get_cumulative_mass(t.edges.right)
            t.edges.set_columns(left=new_left, right=new_right, parent=t.edges.parent, child=t.edges.child)
            new_ts = t.tree_sequence()
            new_ts.dump(ts_path)
        ts_list.append(new_ts)

    # Check individuals are in the same order across tree sequences
    selected_ind = [i for i in ts_list[0].individuals() if int(i.metadata['individual_name']) in proband_pca]
    chr1_order = [i.metadata['individual_name'] for i in selected_ind]
    for chrom in range (2,23):
        selected_ind = [i for i in ts_list[chrom - 1].individuals() if int(i.metadata['individual_name']) in proband_pca]
        chr_order = [i.metadata['individual_name'] for i in selected_ind]
        if chr_order != chr1_order:
            print(f"Warning: individuals in chromosome {chrom} different to chromosome 1.")
    
    # ts rescaled PCA
    mode = 'branch'
    start_time = time.time()
    prin_comp = ts_linop_pca(ts_list, mode=mode)
    end_time = time.time()
    print(end_time - start_time)

    pc_grm_df = pd.DataFrame(prin_comp)
    pc_grm_df['proband'] = [int(i) for i in chr1_order]
    pc_grm_df = pc_grm_df.merge(meta_df, left_on="proband", right_on="ind").drop(columns=["proband"])
    pc_grm_df.to_csv(f"output/balsac/genome/{mode}_pca.csv", index=False)


if __name__ == "__main__":
    main()
