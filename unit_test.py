#!/usr/bin/env python
import numpy as np

# Two input arguments, simulated data and experimental mass data

from execute_process_mass_data import array_pmf_results

standard_dev = []

for line in array_pmf_results:
    standard_dev.append(np.std(line['pmf']))

print(standard_dev)
