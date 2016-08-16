#!/usr/bin/env python

from process_mass_data import *

all_results_import = read_output_data('comp_mass_output_tmp.csv')

filepath = 'settings/process_config.yaml'
data = yaml_loader(filepath)

mof_array = data['mof_array']
mof_densities_import = {}
mof_experimental_mass = {}
for mof in mof_array:
    mof_densities_import.copy()
    mof_densities_import.update({ mof : data['mofs'][mof]['density']})
    mof_experimental_mass.copy()
    mof_experimental_mass.update({ mof : data['mofs'][mof]['experimental_mass'] })
num_mixtures = data['num_mixtures']
stdev = data['stdev']
mrange = data['mrange']
gases = data['gases']

interpolate_data_results = interpolate_data(mof_array, all_results_import, mof_densities_import, gases)
calculate_pmf_results = calculate_pmf(interpolate_data_results, mof_array, mof_experimental_mass, stdev, mrange)
create_bins_results = create_bins(mof_array, calculate_pmf_results, gases)
bin_compositions_results = bin_compositions(gases, mof_array, create_bins_results, calculate_pmf_results)
plot_binned_pmf_array(gases, mof_array, bin_compositions_results, create_bins_results)
