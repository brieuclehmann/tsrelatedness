import numpy as np
import pandas as pd
import time
import tskit
import os

# Download chromosome 21 from https://doi.org/10.5281/zenodo.7702392
# Use tsunzip to unzip the .tsz file

# Subset n individuals
def subset_ts(input_path, n):
    ts = tskit.load(input_path)
    samples = ts.samples()[:(2*n)]
    sub_ts = ts.simplify(samples)
    return sub_ts

def benchmark_GRM(ts):
    before = time.perf_counter()
    D = ts.divergence_matrix(mode="branch")
    duration = time.perf_counter() - before
    data = {"D_time": duration}

    before = time.perf_counter()
    B = ts.genetic_relatedness_matrix(mode="branch")
    duration = time.perf_counter() - before
    data["B_time"] = duration
    return data

def benchmark_GRM_vec(ts):
    before = time.perf_counter()
    w = np.ones(ts.num_samples)
    x = ts.genetic_relatedness_vector(w, mode="branch")
    duration = time.perf_counter() - before
    return {"matvec": duration}


def run_benchmark():
    data = []
    for k in range(1, 7):
        path = f"tmp/chr21_10_{k}.ts"
        if not os.path.exists(path):
            ts = subset_ts("data/simulated_chrom_21.ts", 10**k)
            if not os.path.exists("tmp"):
                os.makedirs("tmp")
            ts.dump(path)
        else:
            ts = tskit.load(path)
        if k == 1:
          print(benchmark_GRM(ts))
        datum = benchmark_GRM_vec(ts)
        datum["n"] = int(ts.num_samples / 2)
        data.append(datum)
        df = pd.DataFrame(data)
        print(df)
        df.to_csv("output/benchmarking/matvec_timing.csv")

if __name__ == "__main__":
    run_benchmark()
