import pandas as pd
import numpy as np

# Define the parameter values
methods = ['eGRM', 'ts.genetic_relatedness_matrix']
modes = ['branch']
sample_sizes = [2**7, 2**8, 2**9, 2**10, 2**11, 2**12]
sequence_lengths = [10**4, 10**5, 10**6, 10**7, 10**8]

# Create a dictionary to hold the parameter values
parameters = {
    'Method': methods,
    'Mode': modes,
    'Number of samples': sample_sizes,
    'Sequence length': sequence_lengths
}

# Set the default values
default_sample_size = 2**10
default_sequence_length = 10**7
mode = 'branch'

# Create a list of dictionaries for each combination of parameter values
data = []
for method in methods:
      for sample_size in sample_sizes:
          data.append({
              'Method': method,
              'Mode': mode,
              'Num_samples': sample_size,
              'Sequence_length': default_sequence_length
          })
      for sequence_length in sequence_lengths:
          data.append({
              'Method': method,
              'Mode': mode,
              'Num_samples': default_sample_size,
              'Sequence_length': sequence_length
          })

# Create the pandas data frame
df = pd.DataFrame(data)
df = df.drop_duplicates()

# Save the data frame
df.to_csv('data/grm_simulation_parameters.csv', index=False)
