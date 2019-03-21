#!/usr/bin/env python
import sys

from datetime import datetime
from mof_array.pmf.process_mass_data import (read_output_data,
                                            yaml_loader,
                                            write_output_data,
                                            import_experimental_results,
                                            import_simulated_data,
                                            calculate_pmf,
                                            create_bins,
                                            bin_compositions,
                                            array_pmf,
                                            plot_binned_pmf_array,
                                            save_array_pmf_data,
                                            save_raw_pmf_data,
                                            information_gain,
                                            choose_best_arrays)

all_results_import = read_output_data(sys.argv[1])
experimental_mass_import = read_output_data(sys.argv[2])

filepath = 'settings/process_config.yaml'
data = yaml_loader(filepath)

mof_array = data['mof_array']
mof_densities_import = {}
mof_experimental_mass = {}

for mof in mof_array:
    mof_densities_import.copy()
    mof_densities_import.update({ mof : data['mofs'][mof]['density']})

num_mixtures = data['num_mixtures']
stdev = data['stdev']
mrange = data['mrange']
gases = data['gases']
number_mofs = data['number_mofs']
number_bins = data['number_bins']

experimental_mass_results, experimental_mass_mofs, experimental_mofs = import_experimental_results(mof_array, experimental_mass_import, mof_densities_import, gases)
import_data_results = import_simulated_data(experimental_mofs, all_results_import, mof_densities_import, gases)
calculate_pmf_results = calculate_pmf(experimental_mass_results, import_data_results, experimental_mofs, stdev, mrange)
array_pmf_results, list_of_arrays = array_pmf(gases, number_mofs, experimental_mofs, calculate_pmf_results, experimental_mass_mofs)
create_bins_results = create_bins(experimental_mofs, calculate_pmf_results, gases, number_bins)
bin_compositions_results = bin_compositions(gases, list_of_arrays, create_bins_results, array_pmf_results, experimental_mass_mofs)
kl_divergence = information_gain(gases, list_of_arrays, bin_compositions_results, create_bins_results)
ordered_by_kld_product, ordered_by_gas, all_arrays_ranked = choose_best_arrays(gases, number_mofs, kl_divergence)

save_raw_pmf_data(calculate_pmf_results)
plot_binned_pmf_array(gases, list_of_arrays, create_bins_results, bin_compositions_results)
save_array_pmf_data(gases, list_of_arrays, create_bins_results, bin_compositions_results)
write_output_data('saved_results/ordered_by_gas_%s.csv' % (datetime.now().strftime("%Y_%m_%d__%H_%M_%S")), ordered_by_gas)
write_output_data('saved_results/ordered_by_kld_product_%s.csv' % (datetime.now().strftime("%Y_%m_%d__%H_%M_%S")), ordered_by_kld_product)
write_output_data('saved_results/all_arrays_ranked_%s.csv' % (datetime.now().strftime("%Y_%m_%d__%H_%M_%S")), all_arrays_ranked)
