import argparse
import sys
import pandas as pd
import tskit

sys.path.insert(0, "utils")
import grm_tools

def main(args):

    ###################
    #### Load data ####
    ###################

    dir_out = 'output/balsac/chr{}/sim{}'.format(
        args.chromosome, args.repetition)

    ts = tskit.load(dir_out + "/trees_recap.ts")

    ######################
    #### Compute GRMs ####
    ######################

    # Load metadata for subsample of 50 individuals from each village
    selected_df = pd.read_csv("data/balsac_subsample.csv")

    # Extract metadata and sample nodes
    grm_ind = [i for i in ts.individuals() if int(i.metadata['individual_name']) in list(selected_df['ind'])]
    col_names = [str(i.metadata['individual_name']) for i in grm_ind]
    for mode in ["branch"]:
        grm_recap = grm_tools.compute_grm(ts, grm_ind, mode = mode)
        grm_recap_df = pd.DataFrame(grm_recap, columns = col_names)
        grm_recap_df.to_csv(dir_out + "/grm_" + mode + "_recap.csv", index=False)
    print("Finished.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-chr", "--chromosome",
        type=int,
        help="specify chromosome number to be simulated"
        )
    parser.add_argument("-rep", "--repetition",
        default=1,
        type=int,
        help="specify repetition"
        )
    args = parser.parse_args()

    main(args)
