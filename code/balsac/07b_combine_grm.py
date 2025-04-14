import argparse
import pandas as pd
import tskit

def main(args):

    ###################
    #### Load data ####
    ###################
    chrom_lengths = pd.read_csv("data/maps/GRCh37_chromosome_length.csv")
    mode = 'branch'

    meta_df = pd.read_csv("genome_simulations/data/balsac_proband_meta.csv")

    chr1_grm_file = f"output/balsac/chr1/sim{args.repetition}/grm_{mode}_recap.csv"
    grm = pd.read_csv(chr1_grm_file)
    col_names = grm.columns.values
    grm_mat = grm.to_numpy() * chrom_lengths.length_bp[0]

    for chrom in range(2,23):
        grm_file = f"output/balsac/chr{chrom}/sim{args.repetition}/grm_{mode}_recap.csv"
        this_grm = pd.read_csv(grm_file)
        this_chrom_length = chrom_lengths.length_bp[chrom - 1]
        grm_mat += this_grm.to_numpy() * this_chrom_length

    ###################
    #### Save data ####
    ###################
    out_file = f"output/balsac/genome/sim{args.repetition}/grm_{mode}_recap.csv"

    full_grm = pd.DataFrame(grm_mat, columns = col_names)
    full_grm.to_csv(out_file, index = False)

    print("Finished.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-rep", "--repetition",
        default=1,
        type=int,
        help="specify repetition"
        )
    args = parser.parse_args()

    main(args)
