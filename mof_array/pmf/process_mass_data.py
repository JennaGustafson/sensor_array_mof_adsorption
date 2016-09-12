"""This code imports mass adsorption data from simulated output,
for multiple MOFs and gas mixtures, and calculates the probability that
specific gases are present in specific mole fractions based on the total mass
adsorbed for each MOF.
"""

from math import isnan, log10
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

def interpolate_data(mofs_list, all_results, mof_densities, gases):
    """Creates additional gas mixtures and calculates probability mass functions (pmf) for all mass values.

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

        all_results_temp = [ row for row in all_results if row['MOF'] == mof]

        for index in range(len(masses)):
            temp_dict = all_results_temp[index].copy()
            temp_dict.update({ 'Mass_mg/cm3' : masses[index] })
            interpolate_results.extend([temp_dict])

    return(interpolate_results)

def add_random_gas(comps, num_mixtures):
    # Adds random gas mixtures to the original data, between min and max of original mole fractions.

    while (len(comps) < 78 + num_mixtures):
        random_gas = ([0.5 * round(random(), 3), 0.5 * round(random(), 3), 0.2 * round(random(), 3)])
        predicted_mass = interp_dat(random_gas)
        if sum(random_gas) <= 1 and not isnan(predicted_mass):
            comps.append(random_gas)
            masses.extend(predicted_mass)

def calculate_pmf(interpolate_data_results, mofs_list, experimental_mass, stdev, mrange):
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

        # Calculates all pmfs based on the experimental mass and normal probability distribution.
        probs = [(ss.norm.cdf(row['Mass_mg/cm3'] + mrange, float(experimental_mass[mof]),
                              stdev * float(experimental_mass[mof])) -
                  ss.norm.cdf(row['Mass_mg/cm3'] - mrange, float(experimental_mass[mof]),
                              stdev * float(experimental_mass[mof]))) for row in
                              interpolate_data_results if row['MOF'] == mof]
        norm_probs = [(i / sum(probs)) for i in probs]

        # Combine mole fractions, mass values, and pmfs into a numpy array for the dictionary creation.
        all_results_temp = [row for row in interpolate_data_results if row['MOF'] == mof]

        for index in range(len(norm_probs)):
            temp_dict = all_results_temp[index].copy()
            temp_dict.update({ 'PMF' : norm_probs[index] })
            pmf_results.extend([temp_dict])

    return(pmf_results)

def create_bins(mofs_list, calculate_pmf_results, gases):
    """Creates bins for all gases, ranging from the lowest to highest mole fractions for each.

    Keyword arguments:
    mofs_list -- list of mofs used in analysis
    calculate_pmf_results -- dictionary output from the calculate_pmf function
    gases -- list of present gases
    """

    # Creates numpy array of all compositions, needed to calculate min/max of each gas's mole frac.
    mof = mofs_list[0]
    temp_one_mof_results = [row for row in calculate_pmf_results if row['MOF'] == mof]
    comps_array = np.array([[float(row[gas]) for gas in gases] for row in temp_one_mof_results])

    bin_range = np.column_stack([np.linspace(min(comps_array[:,index]),
        max(comps_array[:,index]), 12) for index in range(len(gases))])

    bins = [ { gases[index] : row[index] for index in range(len(gases)) } for row in bin_range]

    return(bins)

def bin_compositions(gases, mof_array, create_bins_results, calculate_pmf_results):
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
        for gas_name in gases:
            binned_data =[]
            binned_probability_temporary = []

            # Assigns pmf to bin value (dictionary) by checking whether mole frac is
            # between the current and next bin value.
            for row in calculate_pmf_results:
                 for i in range(1, len(create_bins_results)):
                    if ( float(row[gas_name]) >= create_bins_results[i - 1][gas_name] and
                         float(row[gas_name]) < create_bins_results[i][gas_name] and
                         row['MOF'] == mof_name
                       ):
                        binned_data.append({ 'probability': row['PMF'],
                                             'bin': create_bins_results[i - 1][gas_name]
                                           })

            # Loops through all of the bins and averages the pmfs into their assgned bins.
            for b in create_bins_results:
                average = []
                for line in binned_data:
                     if b[gas_name] == line['bin']:
                        average.append(line['probability'])
                if average == []:
                    binned_probability_temporary.append({'bin' : line['bin'],
                        'average probability' : 0})
                else:
                    binned_probability_temporary.append({'bin' : line['bin'],
                        'average probability' : np.mean(average)})

            # Creates list of binned probabilities in order to loop through and normalize, sum must be 1.
            temporary_pmf_list = [row['average probability'] for row in binned_probability_temporary]
            normalized_temporary_pmf = [number / sum(temporary_pmf_list) for number in
                temporary_pmf_list]
            binned_probability.extend([{'mof' : mof_name, 'gas' : gas_name,
                'bin' : binned_probability_temporary[i]['bin'],'average probability' :
                normalized_temporary_pmf[i]} for i in range(0, len(normalized_temporary_pmf))])
    return(binned_probability)

def normalize_binned_pmf(gas_names, number_mofs, mof_names, bin_compositions_results, create_bins_results):
    """Normalized the binned probability mass functions for a MOF array"""
    # number_mofs will conatin a lower limit # mofs and an upper limit
    num_mofs = min(number_mofs)
    mof_array_list = []
    while num_mofs < max(number_mofs):
        mof_array_list.extend(list(combinations(mof_names, num_mofs)))
        num_mofs += 1

    normalized_pmf = []
    for mof_array_temp in mof_array_list:
        for gas_name in gas_names:
            compound_pmfs = []
            for mof in mof_array_temp: # set up names for each mof array based on # mofs and combos
                if compound_pmfs == []:
                    compound_pmfs = np.array([point['average probability'] for point in
                        bin_compositions_results if point['mof'] == mof and point['gas'] == gas_name])
                else:
                    compound_pmfs *= np.array([point['average probability'] for point in
                        bin_compositions_results if point['mof'] == mof and point['gas'] == gas_name])
            normalized_compound_pmfs = [number / sum(compound_pmfs) for number in compound_pmfs]
            normalized_pmf.append({'mof array': tuple(mof_array_temp), 'gas' : gas_name, 'pmf' :
             normalized_compound_pmfs})
    return(normalized_pmf)

def plot_binned_pmf_array(gas_names, mof_names, bin_compositions_results, create_bins_results):
    """Calculates compound pmfs for MOF array and plots vs mole fraction for each gas.

    Keyword arguments:
    gas_names -- list of gases specified by user
    mof_names -- list of MOFs in array, specified by user
    bin_compositions_results -- dictionary result from bin_compositions function
    create_bins_results -- dictionary result from create_bins
    """

    # plot_PMF = plt.figure()
    # plt.plot([b[gas_name] for b in create_bins_results], [point for point in
    #     normalized_compound_pmfs], 'bo')
    # plt.savefig("%s_plot_PMF_%s_%s.png" % (datetime.now().strftime("%Y_%m_%d__%H_%M_%S"), str(gas_name) , "_".join(mof_names)))
    # plt.close(plot_PMF)

def information_gain(normalize_binned_pmf_results, create_bins_results):
    """Calculates the Kullback-Liebler Divergence of a MOF array with each gas component.

    Keyword arguments:
    normalize_binned_pmf_results -- list of dictionaries including array names, gases, pmfs
    create_bins_results -- dictionary result from create_bins
    """

    array_gas_info_gain = []
    reference_prob = 1/len(create_bins_results)

    for array_pmf in normalize_binned_pmf_results:
        kl_divergence = sum([float(pmf)*log10(float(pmf)/reference_prob) for pmf in array_pmf['pmf'] if pmf != 0])
        array_gas_info_gain.append({'mof array': array_pmf['mof array'], 'gas': array_pmf['gas'], 'KLD': round(kl_divergence,4)})

    return(array_gas_info_gain)
