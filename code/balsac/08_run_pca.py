import os
import argparse
import sys
import time

import pandas as pd
import numpy as np
import tskit
import stdpopsim    
print(stdpopsim.__version__)

from scipy.linalg import eigh
import matplotlib
import matplotlib.pyplot as plt

sys.path.insert(0, "utils")
import pedigree_tools
import grm_tools

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
    ### Pedigree PCA ###
    ####################
    n_pc = 8
    prm_df = pd.read_csv("output/chr3/prm_full.csv")

    prm_ind = prm_df.columns.values
    prm_df = prm_df.set_index(prm_ind)
    region_df = region_df.astype({"ind":"str"})
    region_map = pd.Series(region_df.proband_region.values,index=region_df.ind).to_dict()
    prm_regions = [region_map[i] for i in prm_ind]

    # Filter out probands with shallow average founder depth
    prm_df = prm_df.loc[[str(i) for i in proband_pca], [str(i) for i in proband_pca]]
    prm_ind = prm_df.columns.values

    prm_mat = prm_df.to_numpy()
    n_ind = prm_mat.shape[0]
    prm_eigval, prm_eigvec = eigh(prm_mat, subset_by_index=[n_ind-n_pc-1, n_ind-1])
    pc_prm = (prm_mat @ prm_eigvec[:, ::-1]) 
    pc_prm_df = pd.DataFrame(pc_prm)
    pc_prm_df['proband'] = [int(i) for i in prm_ind]
    pc_prm_df = pc_prm_df.merge(meta_df, left_on = "proband", right_on="ind").drop(columns=["proband"])
    pc_prm_df.to_csv("output/balsac/pedigree_pca.csv", index=False)

    ####################
    #### Branch PCA ####
    ####################
    for sim in [1]:
        ts_dir_out = f'output/balsac/chr3/sim{sim}'

        ts = tskit.load(ts_dir_out + "/trees_recap.ts")

        ### Rescale tree sequence ###
        ts_path = ts_dir_out + "/trees_recap_rescaled.ts"
        if os.path.exists(ts_path):
            new_ts = tskit.load(ts_path)
        else:
            species = stdpopsim.get_species("HomSap")
            contig = species.get_contig("chr3", genetic_map="HapMapII_GRCh37")
            rate_map = contig.recombination_map
            t = ts.dump_tables()
            t.mutations.clear()
            t.sites.clear()
            new_left = rate_map.get_cumulative_mass(t.edges.left)
            new_right = rate_map.get_cumulative_mass(t.edges.right)
            t.edges.set_columns(left=new_left, right=new_right, parent=t.edges.parent, child=t.edges.child)
            new_ts = t.tree_sequence()
            new_ts.dump(ts_path)

        selected_ind = [i for i in new_ts.individuals() if int(i.metadata['individual_name']) in proband_pca]
        pca_ind = [i.id for i in selected_ind]
        balsac_ind = [int(i.metadata['individual_name']) for i in selected_ind]
    
        for mode in ['branch']:
        # ts rescaled PCA
            print("Number of PC individuals :", len(pca_ind))
            start_time = time.time()
            out = new_ts.pca(individuals=np.asarray(pca_ind), num_iterations=5, random_seed=1, num_components=8)
            end_time = time.time()
            print(end_time - start_time)
            #grm_trace = np.trace(grm_rescaled)
            prin_comp = out.factors
            pc_grm_df = pd.DataFrame(prin_comp)
            pc_grm_df['proband'] = [int(i) for i in balsac_ind]
            pc_grm_df = pc_grm_df.merge(meta_df, left_on="proband", right_on="ind").drop(columns=["proband"])
            pc_grm_df.to_csv(f"output/balsac/chr3/sim{sim}/{mode}_pca.csv", index=False)

if __name__ == "__main__":
    main()
