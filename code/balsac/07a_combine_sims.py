import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

prm_file = "output/chr3_prm.csv"
prm_df = pd.read_csv(prm_file)

prm_df['proband1'] = prm_df.columns.values
prm_long = pd.melt(prm_df, id_vars='proband1', var_name='proband2', value_name='prm')
prm_long = prm_long[prm_long['proband1'] <= prm_long['proband2']]

# Get selected individuals in region order
selected_df = pd.read_csv("data/balsac_subsample.csv")
order_cols = selected_df.sort_values("proband_region")['ind'].tolist()
order_cols = [str(i) for i in order_cols]

# Get relative information
relative_df = pd.read_csv("data/balsac_relatives.csv")
relative_df = relative_df.loc[relative_df['proband1'].isin(selected_df['ind']) & relative_df['proband2'].isin(selected_df['ind'])]
relative_df = relative_df[['proband1', 'proband2', 'generation', 'relationship']]

grm_dir = "grms/balsac/chr3/"
full_df = pd.DataFrame()
for type in ["branch_recap"]:
    for sim in range(1,101):
        grm_file = grm_dir + "sim" + str(sim) + "/" + type + ".csv"
        try:
            grm_df = pd.read_csv(grm_file)
            grm_df['proband1'] = grm_df.columns.values
            # Reset index to make row names a column
            melted = pd.melt(grm_df, id_vars='proband1', var_name='proband2', value_name='grm')
            melted = melted[melted['proband1'] <= melted['proband2']]
            melted['sim'] = sim
            melted['type'] = type
            full_df = pd.concat([full_df, melted])
        except:
            pass

full_df = full_df.merge(prm_long, on = ["proband1", "proband2"])
full_df = full_df.astype({"proband1":"int64", "proband2":"int64"})
full_df = full_df.merge(relative_df, on = ["proband1", "proband2"], how = "left")

full_df.to_csv("output/relatives_relatedness.csv", index=False)


