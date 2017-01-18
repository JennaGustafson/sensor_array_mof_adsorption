"""This code imports mass adsorption data and "experimental" adsorption data
from simulated output, for multiple MOFs and gas mixtures, and calculates the
probability that specific gases are present in specific mole fractions
(in each experimental case) based on the total mass adsorbed for each MOF.
"""

from math import isnan, log
import csv
from random import random
from itertools import combinations

import yaml
import numpy as np
from matplotlib import pyplot as plt
import scipy.stats as ss
from scipy.spatial import Delaunay
import scipy.interpolate as si
from datetime import datetime

# Function imports csv file as a dictionary

def read_output_data(filename):
    with open(filename,newline='') as csvfile:
        output_data = csv.DictReader(csvfile, delimiter="\t")
        return list(output_data)

def yaml_loader(filepath):
    with open(filepath, 'r') as yaml_file:
        data = yaml.load(yaml_file)
    return(data)

def import_experimental_results(mofs_list, experimental_mass_import, mof_densities, gases):
    """Imports the experimental data and puts it in dictionary format

    Keyword arguments:
    mofs_list -- list of MOF structures simulated
    experimental_mass_import -- dictionary formatted experimental results for each mof
    mof_densities -- dictionary of densities
    gases -- list of gases in simulated mixtures
    """
    experimental_results = []
    experimental_mass_mofs = []

    for mof in mofs_list:

        # Calculates masses in terms of mg/(cm3 of framework)
        masses = [float(mof_densities[row['MOF']]) * float(row['Mass']) for row in
                    experimental_mass_import if row['MOF'] == mof]

        # Saves composition values as a list, necessary type for the Delaunay input argument
        comps = []
        for row in experimental_mass_import:
            if row['MOF'] == mof:
                comps.extend([[float(row[gas]) for gas in gases]])

        # List of experimental masses for current mof
        experimental_mass_temp = [ row for row in experimental_mass_import if row['MOF'] == mof]

        for index in range(len(masses)):
            temp_dict = experimental_mass_temp[index].copy()
            temp_dict.update({ 'Mass_mg/cm3' : round(masses[index], 6) })
            experimental_results.extend([temp_dict])

        # Dictionary format of all the experimental data
        temp_list = {'MOF' : mof, 'Mass' :[row['Mass_mg/cm3'] for row in experimental_results if row['MOF'] == mof]}
        experimental_mass_mofs.append(temp_list)

    return(experimental_results, experimental_mass_mofs)

def import_interpolate_data(mofs_list, all_results, mof_densities, gases):
    """Imports simulated data and puts it in dictionary format
    If desired, interpolation is performed and may be used to create a denser data set

    Keyword arguments:
    mofs_list -- names of all mofs
    all_results -- imported csv of outputs as dictionary
    mof_densities -- dictionary of densities for each mof
    gases -- list of all gases
    """
    interpolate_results = []
    for mof in mofs_list:
        # Calculates masses in terms of mg/(cm3 of framework)
        masses = [float(mof_densities[row['MOF']]) * float(row['Mass']) for row in
                    all_results if row['MOF'] == mof]

        # Saves composition values as a list, necessary type for the Delaunay input argument
        comps = []
        for row in all_results:
            if row['MOF'] == mof:
                comps.extend([[float(row[gas]) for gas in gases]])
        d = Delaunay(np.array(comps)[:,range(len(gases)-1)])
        interp_dat = si.LinearNDInterpolator(d, masses)

        # Update dictionary with simulated data in terms of mg/cm3
        all_results_temp = [ row for row in all_results if row['MOF'] == mof]
        for index in range(len(masses)):
            temp_dict = all_results_temp[index].copy()
            temp_dict.update({ 'Mass_mg/cm3' : round(masses[index], 6) })
            interpolate_results.extend([temp_dict])

    return(interpolate_results)

def add_random_gas(comps, num_mixtures):
    """ Adds random gas mixtures to the original data, between min and max of original mole fractions

    Keyword arguments:
    comps -- all simulated gas compositions
    num_mixtures -- specify integer number of mixtures to add
    """
    while (len(comps) < 78 + num_mixtures):
        random_gas = ([0.5 * round(random(), 3), 0.5 * round(random(), 3), 0.2 * round(random(), 3)])
        predicted_mass = interp_dat(random_gas)
        if sum(random_gas) <= 1 and not isnan(predicted_mass):
            comps.append(random_gas)
            masses.extend(predicted_mass)

def calculate_pmf(experimental_mass_results, interpolate_data_results, mofs_list, experimental_mass, stdev, mrange):
    """Calculates probability mass function of each data point

    Keyword arguments:
    interpolate_data_results -- dictionary, results of interpolate_data method
    mofs_list -- names of all MOFs
    experimental_mass -- masses "experimentally" obtained for each MOF
    stdev -- standard deviation for the normal distribution
    mrange -- range for which the difference between cdfs is calculated
    """
    pmf_results = []
    for mof in mofs_list:

        mof_temp_dict = []
        # Combine mole fractions, mass values, and pmfs into a numpy array for the dictionary creation.
        all_results_temp = [row for row in interpolate_data_results if row['MOF'] == mof]

        #Loop through all of the experimental masses for each MOF, read in and save comps
        experimental_mass_data = [data_row['Mass_mg/cm3'] for data_row in experimental_mass_results
                                    if data_row['MOF'] == mof]
        # for experimental_mass in experimental_mass_results if experimental_mass['MOF'] == mof:
        # mof_mass = experimenatal_mass['Mass_mg/cm3']
        for mof_mass in experimental_mass_data:

            new_temp_dict = []
            # Calculates all pmfs based on the experimental mass and normal probability distribution.
            probs = [(ss.norm.cdf(row['Mass_mg/cm3'] + mrange, float(mof_mass),
                              stdev * float(mof_mass)) -
                      ss.norm.cdf(row['Mass_mg/cm3'] - mrange, float(mof_mass),
                              stdev * float(mof_mass))) for row in
                              interpolate_data_results if row['MOF'] == mof]
            norm_probs = [(i / sum(probs)) for i in probs]

            if mof_temp_dict == []:
                for index in range(len(norm_probs)):
                    mof_temp_dict = all_results_temp[index].copy()
                    mof_temp_dict.update({ 'PMF_%s' % str(round(mof_mass, 2)) : norm_probs[index] })
                    new_temp_dict.extend([mof_temp_dict])
                mass_temp_dict = new_temp_dict
            else:
                for index in range(len(norm_probs)):
                    mof_temp_dict = mass_temp_dict[index].copy()
                    mof_temp_dict.update({ 'PMF_%s' % str(round(mof_mass, 2)) : norm_probs[index] })
                    new_temp_dict.extend([mof_temp_dict])
                mass_temp_dict = new_temp_dict

        pmf_results.extend(mass_temp_dict)

    return(pmf_results)

def create_bins(mofs_list, calculate_pmf_results, gases):
    """Creates bins for all gases, ranging from the lowest to highest mole fractions for each.

    Keyword arguments:
    mofs_list -- list of mofs used in analysis
    calculate_pmf_results -- dictionary output from the calculate_pmf function
    gases -- list of present gases
    """
    num_bins = 12

    # Creates numpy array of all compositions, needed to calculate min/max of each gas's mole frac.
    mof = mofs_list[0]
    temp_one_mof_results = [row for row in calculate_pmf_results if row['MOF'] == mof]
    comps_array = np.array([[float(row[gas]) for gas in gases] for row in temp_one_mof_results])

    bin_range = np.column_stack([np.linspace(min(comps_array[:,index]), max(comps_array[:,index]) +
        (max(comps_array[:,index])-min(comps_array[:,index]))/num_bins, num_bins + 1) for index in range(len(gases))])

    bins = [ { gases[index] : row[index] for index in range(len(gases)) } for row in bin_range]

    return(bins)

def bin_compositions(gases, mof_array, create_bins_results, calculate_pmf_results, experimental_mass_mofs):
    """Sorts pmfs into bins created by create_bins function.

    Keyword arguments:
    gases -- list of gases specified as user input
    mof_array -- list of MOFs in the array, specified as user input
    create_bins_results -- dictionary containing bins for each gas
    interpolate_pmf_results -- dictionary of pmfs associated with MOFs/gases
    as a result of interpolate_pmf function
    """

    binned_probability = []
    for mof_name in mof_array:
        experimental_mass = [row['Mass'] for row in experimental_mass_mofs if row['MOF'] == mof_name]
        for gas_name in gases:
            binned_data = []

            # Assigns pmf to bin value (dictionary) by checking whether mole frac is
            # between the current and next bin value.
            for row in calculate_pmf_results:
                 for i in range(1, len(create_bins_results)):
                    if ( float(row[gas_name]) >= create_bins_results[i - 1][gas_name] and
                         float(row[gas_name]) < create_bins_results[i][gas_name] and
                         row['MOF'] == mof_name
                       ):
                        for mof_mass in experimental_mass[0]:
                            binned_data.append({ 'probability_%s' % str(round(mof_mass,2)): row['PMF_%s' % str(round(mof_mass,2))],
                                                 'bin': create_bins_results[i - 1][gas_name]
                                               })

            # Loops through all of the bins and averages the pmfs into their assgned bins.
            for mof_mass in experimental_mass[0]:
                binned_probability_temporary = []
                for b in create_bins_results[0:len(create_bins_results)-1]:
                    average = []
                    for line in binned_data:
                         if b[gas_name] == line['bin'] and 'probability_%s' % str(round(mof_mass,2)) in line.keys():
                            average.append(line[ 'probability_%s' % str(round(mof_mass,2))])
                    if average == []:
                        binned_probability_temporary.append({'bin' : b[gas_name],
                            'average probability' : 0})
                    else:
                        binned_probability_temporary.append({'bin' : b[gas_name],
                            'average probability' : np.mean(average)})

                # Creates list of binned probabilities in order to loop through and normalize, sum must be 1.
                temporary_pmf_list = [row['average probability'] for row in binned_probability_temporary]
                normalized_temporary_pmf = [number / sum(temporary_pmf_list) for number in
                    temporary_pmf_list]
                binned_probability.extend([{'mof' : mof_name, 'gas' : gas_name,
                    'bin' : binned_probability_temporary[i]['bin'],'average probability_%s' % str(round(mof_mass,2)) :
                    normalized_temporary_pmf[i]} for i in range(0, len(normalized_temporary_pmf))])

    return(binned_probability)


def compound_pmf_for_mof_array(mof_array, experimental_mass_mofs, gas_name, mof_mass_index, bin_compositions_results):
    """Combines pmfs for a mof array and gas combination, used in method 'normalize_binned_pmf'

    Keyword arguments:
    mof_array -- list of mofs in array
    experimental_mass_mofs -- ordered list of dictionaries with each experimental mof/mass
    gas_name -- current gas
    mof_mass_index -- experimental mass integer identifier
    bin_compositions_results -- list of dictionaries including mof, gas, bin, averaged probability
    """
    compound_pmfs = None
    for mof in mof_array:
        # Call the experimental mass for one mof and experiment, labeled index here
        # Take [0] to retrieve the mass as a float (from a list with one element)
        mof_mass = [ row['Mass'][mof_mass_index] for row in experimental_mass_mofs if row['MOF'] == mof ][0]
        key = 'average probability_%s' % str(round(mof_mass,2))
        mof_gas_pmf = []
        for point in bin_compositions_results:
            # Creates list of pmf values for a MOF/gas/exp combination
            if point['mof'] == mof and point['gas'] == gas_name and key in point.keys():
                mof_gas_pmf.append(point['average probability_%s' % str(round(mof_mass,2))])
        # Save list as numpy array for joint prob calculation
        mof_gas_pmf = np.array(mof_gas_pmf)

        # Joint prob, multiply pmf values elementwise
        if compound_pmfs is not None:
            compound_pmfs *= mof_gas_pmf
        else:
            compound_pmfs = mof_gas_pmf
    # Normalize joint probability, sum of all points is 1
    normalized_compound_pmfs = [ number / sum(compound_pmfs) for number in compound_pmfs ]
    return normalized_compound_pmfs

def normalize_binned_pmf(gas_names, number_mofs, mof_names, bin_compositions_results, experimental_mass_mofs):
    """Normalizes the binned probability mass functions for a MOF array

    Keyword arguments:
    gas_names -- list of gases
    number_mofs -- lower and upper limit of desired number of mofs in array
    mof_names -- list of all mofs
    bin_compositions_results -- list of dictionaries including mof, gas, bin, averaged probability
    experimental_mass_mofs -- ordered list of dictionaries with each experimental mof/mass
    """
    num_mofs = min(number_mofs)
    mof_array_list = []
    # Creates list of MOF arrays, all combinations from min to max number of MOFs
    while num_mofs <= max(number_mofs):
        mof_array_list.extend(list(combinations(mof_names, num_mofs)))
        num_mofs += 1

    normalized_pmf = []
    # Nested loops take all combinations of array/gas/experiment
    for mof_array in mof_array_list:
        for gas_name in gas_names:
            for mof_mass_index in range(0,len(experimental_mass_mofs[0]['Mass'])):
                # Calls outside function to calculate joint probability
                normalized_compound_pmfs = compound_pmf_for_mof_array(mof_array, experimental_mass_mofs,
                                                                      gas_name, mof_mass_index,
                                                                      bin_compositions_results
                                                                      )
                normalized_pmf.append({
                    'mof array': tuple(mof_array),
                    'gas' : gas_name,
                    'pmf_%s' % mof_mass_index : normalized_compound_pmfs
                })

    return(normalized_pmf)

def plot_binned_pmf_array(gas_names, mof_names, create_bins_results, experimental_mass_mofs, normalize_binned_pmf_results):
    """Plots pmf vs mole fraction for each gas/MOF array combination

    Keyword arguments:
    gas_names -- list of gases specified by user
    mof_names -- list of MOFs in array, specified by user
    bin_compositions_results -- dictionary result from bin_compositions function
    create_bins_results -- dictionary result from create_bins
    """

    # Identifies experiment of interest
    mof_mass_index = 0
    # All results for specified experiment collected in list of dictionaries
    pmf_per_experiment = [pmf for pmf in normalize_binned_pmf_results if 'pmf_%s' % mof_mass_index in pmf.keys()]
    for each_array_gas_combo in pmf_per_experiment:
        # Y-axis, list of pmf values to plot
        pmfs_to_plot = each_array_gas_combo['pmf_%s' % mof_mass_index]
        gas_name = each_array_gas_combo['gas']
        mof_names = each_array_gas_combo['mof array']
        # X-axis, list of mole fracs to plot, for relevant gas
        comps_to_plot = [b[gas_name] for b in create_bins_results][:len(create_bins_results)-1]
        # Plot and save figure in a directory 'figures'
        plot_PMF = plt.figure()
        plt.plot(comps_to_plot, pmfs_to_plot, 'bo')
        plt.title('Experiment %s, Array %s, Gas %s' % (mof_mass_index, "_".join(mof_names), str(gas_name)))
        plt.savefig("figures/%s_%s_%s_%s.png" % (mof_mass_index, "_".join(mof_names), str(gas_name), datetime.now().strftime("%Y_%m_%d__%H_%M_%S")))
        plt.close(plot_PMF)

def information_gain(normalize_binned_pmf_results, create_bins_results, experimental_mass_mofs):
    """Calculates the Kullback-Liebler Divergence of a MOF array with each gas component.

    Keyword arguments:
    normalize_binned_pmf_results -- list of dictionaries including array names, gases, pmfs
    create_bins_results -- dictionary result from create_bins
    experimental_mass_mofs -- ordered list of dictionaries with each experimental mof/mass
    """
    array_gas_info_gain = []
    reference_prob = 1/len(create_bins_results)
    for mof_mass_index in range(0,len(experimental_mass_mofs[0]['Mass'])):
        # For each experiment, take list of dictionaries with results
        pmf_per_experiment = [pmf for pmf in normalize_binned_pmf_results if 'pmf_%s' % mof_mass_index in pmf.keys()]
        for array_pmf in pmf_per_experiment:
            # For each array/gas combination, calculate the kld
            kl_divergence = sum([float(pmf)*log(float(pmf)/reference_prob,2) for pmf in array_pmf['pmf_%s' % mof_mass_index] if pmf != 0])
            # Result is list of dicts, experiment expressed as integer (1,2,3...), dropping the pmf values
            array_gas_info_gain.append({'experiment': mof_mass_index +1, 'mof array': array_pmf['mof array'],
                                        'gas': array_pmf['gas'], 'KLD': round(kl_divergence,4)})

    return(array_gas_info_gain)

def choose_best_arrays(gas_names, information_gain_results):
    """Choose the best MOF arrays by selecting the top KL scores for each gas

    Keyword arguments:
    gas_names -- list of gases
    information_gain_results -- list of dictionaries including each experiment, mof array, gas, and corresponding kld
    """
    # Sort results from highest to lowest KLD values
    ordered_by_kld = sorted(information_gain_results, key=lambda k: k['KLD'], reverse=True)
    best_overall_array = ordered_by_kld[0]

    # Number of array/experiment combinations for each gas
    num_points_per_gas = int(len(information_gain_results)/len(gas_names))
    best_gas_index = 0
    best_per_gas = []

    # Already sorted from highest to lowest KLD, sort by gas
    best_by_gas = sorted(ordered_by_kld, key=lambda k: k['gas'])

    # Take the highest KLD values for every gas, increasing the loop by number of results per gas
    while best_gas_index < len(information_gain_results):
        best_per_gas.append(best_by_gas[best_gas_index])
        best_gas_index += num_points_per_gas

    # Prints the best arrays for each gas, regardless of experiment number
    for result in best_per_gas:
        print("The best array for %s consists of MOFs: %s" % (result['gas'], result['mof array']))

    # Print all arrays in order of KLD values
    for each in ordered_by_kld:
        print(each)
