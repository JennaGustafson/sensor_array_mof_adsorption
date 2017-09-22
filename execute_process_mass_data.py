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
calculate_pmf_results = calculate_pmf(experimental_mass_results, import_data_results, experimental_mofs, mof_experimental_mass, stdev, mrange)
# save_raw_pmf_data(calculate_pmf_results)

array_pmf_results, labeled_exp_mass_mofs, list_of_arrays = array_pmf(gases, number_mofs, experimental_mofs, calculate_pmf_results, experimental_mass_mofs)
create_bins_results = create_bins(experimental_mofs, calculate_pmf_results, gases, number_bins)
bin_compositions_results = bin_compositions(gases, list_of_arrays, create_bins_results, array_pmf_results, experimental_mass_mofs)
# plot_binned_pmf_array(gases, list_of_arrays, create_bins_results, bin_compositions_results)
# save_array_pmf_data(gases, list_of_arrays, create_bins_results, bin_compositions_results)
kl_divergence = information_gain(gases, list_of_arrays, bin_compositions_results, create_bins_results)
# combined_kld, ordered_by_gas, best_arrays, ordered_kld_w_array, average_kld = choose_best_arrays(gases, kl_divergence)

# write_output_data('saved_results/ordered_by_kld_product_%s.csv' % (datetime.now().strftime("%Y_%m_%d__%H_%M_%S")), combined_kld)
# write_output_data('saved_results/ordered_by_gas_kld_%s.csv' % (datetime.now().strftime("%Y_%m_%d__%H_%M_%S")), ordered_by_gas)
#
# # Print results, including the "best" MOF array structures
# print(' ================ RESULTS ===============')
# print("\nThe MOF array with the highest informaion content, of %s, is: %s \n" % (str(ordered_kld_w_array[0]['KLD']),
#                                                                                  str(ordered_kld_w_array[0]['mof array'])))
# print("\nThe MOF array with the lowest informaion content, of %s, is: %s \n" % (str(ordered_kld_w_array[len(ordered_kld_w_array)-1]['KLD']),
#                                                                                 str(ordered_kld_w_array[len(ordered_kld_w_array)-1]['mof array'])))
#
# for result in best_arrays:
#     print("The best array for %s consists of MOFs: %s" % (result['gas'], result['mof array']))
# print('\nThe average Kullback-Liebler Divergence for a: ')
# for num_mofs in range(0,len(average_kld)):
#     print('\t%s MOF array is %s' %(num_mofs + 1, average_kld[num_mofs]))
