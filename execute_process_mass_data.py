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

gases = ['N2', 'CH4', 'CO2', 'C2H6']
mof_array = ['IRMOF-1', 'HKUST-1', 'NU-125']
interpolate_pmf_results = interpolate_pmf(mofs_import, all_results_import, mof_experimental_mass, mof_densities_import, num_mixtures, stdev, mrange)
create_bins_results = create_bins(interpolate_pmf_results)
bin_compositions_results = bin_compositions(gases, mof_array, create_bins_results, interpolate_pmf_results)
plot_binned_pmf_array(gases, mof_array, bin_compositions_results, create_bins_results)
