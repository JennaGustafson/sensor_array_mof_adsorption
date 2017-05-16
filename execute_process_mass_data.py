#!/usr/bin/env python
import sys

from mof_array.pmf.process_mass_data import (read_output_data,
                                            yaml_loader,
                                            import_experimental_results,
                                            import_simulated_data,
                                            calculate_pmf,
                                            create_bins,
                                            bin_compositions,
                                            array_pmf,
                                            plot_binned_pmf_array,
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
create_bins_results = create_bins(experimental_mofs, calculate_pmf_results, gases, number_bins)
bin_compositions_results = bin_compositions(gases, experimental_mofs, create_bins_results, calculate_pmf_results, experimental_mass_mofs)
array_pmf_results, labeled_exp_mass_mofs = array_pmf(gases, number_mofs, experimental_mofs, bin_compositions_results, experimental_mass_mofs)
plot_binned_pmf_array(gases, experimental_mofs, create_bins_results, array_pmf_results)
kl_divergence = information_gain(array_pmf_results, create_bins_results, labeled_exp_mass_mofs)
best_arrays, ordered_kld_w_array, average_kld = choose_best_arrays(gases, kl_divergence)

# Print results, including the "best" MOF array structures
print(' ================ RESULTS ===============')
print("\nThe MOF array with the highest informaion content, of %s, is: %s \n" % (str(ordered_kld_w_array[0]['KLD']),
                                                                                 str(ordered_kld_w_array[0]['mof array'])))
print("\nThe MOF array with the lowest informaion content, of %s, is: %s \n" % (str(ordered_kld_w_array[len(ordered_kld_w_array)-1]['KLD']),
                                                                                str(ordered_kld_w_array[len(ordered_kld_w_array)-1]['mof array'])))

for result in best_arrays:
    print("The best array for %s consists of MOFs: %s" % (result['gas'], result['mof array']))
print('\nThe average Kullback-Liebler Divergence for a: ')
for num_mofs in range(0,len(average_kld)):
    print('\t%s MOF array is %s' %(num_mofs + 1, average_kld[num_mofs]))
