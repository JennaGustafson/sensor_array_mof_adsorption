#!/usr/bin/env python

from process_mass_data import *

all_results_import = read_output_data('comp_mass_output_tmp.csv')

filepath = 'process_config.yaml'
data = yaml_loader(filepath)

mof_densities_import = data['density']
mofs_import = data['mofs']
mof_experimental_mass = data['experimental_mass']
num_mixtures = data['num_mixtures']
stdev = data['stdev']
mrange = data['mrange']
gases = data['gases']
mof_array = data['mof_array']

interpolate_data_results = interpolate_data(mofs_import, all_results_import, mof_densities_import, gases)
calculate_pmf_results = (interpolate_data_results, mofs_import, mof_experimental_mass, stdev, mrange)
create_bins_results = create_bins(calculate_pmf_results)
bin_compositions_results = bin_compositions(gases, mof_array, create_bins_results,
    interpolate_pmf_results)
plot_binned_pmf_array(gases, mof_array, bin_compositions_results, create_bins_results)
