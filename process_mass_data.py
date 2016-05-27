"""This code imports mass adsorption data from simulated output,
for multiple MOFs and gas mixtures, and calculates the probability that
specific gases are present in specific mole fractions based on the total mass
adsorbed for each MOF.
"""

from math import isnan
import csv
from random import random

import numpy as np
from matplotlib import pyplot as plt
import scipy.stats as ss
from scipy.spatial import Delaunay
import scipy.interpolate as si

# Function imports a tab delimited csv file as dictionary

def mof_density(filename):
    with open(filename, newline='') as csvfile:
        density = csv.DictReader(csvfile, delimiter="\t")
        return list(density)[0]

# Function imports csv file as a dictionary

def read_output_data(filename):
    with open(filename,newline='') as csvfile:
        output_data = csv.DictReader(csvfile, delimiter="\t")
        return list(output_data)

# Reads in csvfile as a list

def mof_names(filename):
    with open(filename,newline='') as csvf:
        names = csv.reader(csvf,delimiter='\n')
        return list(names)


# num_mixtures = 10 #specify how many *additional* mixtures one wants in the analysis (to be interpolated)
# r = .05           #the range for which the pmf is calculated, (cdf mass+r minus cdf mass-r)
# stdev = 0.1       #standard deviation of the  experimental masses for normal distribution curves

def interpolate_pmf(mofs_list, all_results, experimental_mass, mof_densities):
    """Creates additional gas mixtures and calculates probability mass functions (pmf) for all mass values.

    Keyword arguments:
    mofs_list -- names of all mofs
    all_results -- imported csv of outputs as dictionary
    experimental_mass -- dictionary of masses for each mof
    mof_densities -- dictionary of densities for each mof
    """

    pmf_results = []
    for mof in mofs_list:

        # Calculates masses in terms of mg/(cm3 of framework)
        masses = [float(mof_densities[row['MOF']]) * float(row['Mass 1bar']) for row in all_results if row['MOF'] == mof]

        # Saves composition values as a list, necessary type for the Delaunay input argument
        comps = [[float(row['CH4']), float(row['CO2']), float(row['C2H6'])] for row in all_results if row['MOF'] == mof]
        d = Delaunay(comps)
        interp_dat = si.LinearNDInterpolator(d, masses)

        # Adds random gas mixtures to the original data, between min and max of original mole fractions.
        while (len(comps) < 78 + num_mixtures):
            random_gas = ([0.5 * round(random(), 3), 0.5 * round(random(), 3), 0.2 * round(random(), 3)])
            predicted_mass = interp_dat(random_gas)
            if sum(random_gas) <= 1 and not isnan(predicted_mass):
                comps.append(random_gas)
                masses.extend(predicted_mass)

        # Calculates all pmfs based on the experimental mass and normal probability distribution.
        probs = [(ss.norm.cdf(mass + r, float(experimental_mass[mof]),
                              stdev * float(experimental_mass[mof])) -
                  ss.norm.cdf(mass - r, float(experimental_mass[mof]),
                              stdev * float(experimental_mass[mof]))) for mass in masses]
        norm_probs = [(i / sum(probs)) for i in probs]

        # Combine mole fractions, mass values, and pmfs into a numpy array for the dictionary creation.
        comps_mass_prob = np.column_stack((comps, masses, norm_probs))
        pmf_results.extend([{'mof': mof, 'CH4': row[0], 'CO2': row[1],
                             'C2H6': row[2], 'N2' : 1 - (row[0] + row[1] + row[2]),
                             'Mass 1bar' : row[3], 'PMF 1bar' : row[4]} for row in comps_mass_prob])
    return(pmf_results)

def create_bins(interpolate_pmf_results):

    """Creates bins for all gases, ranging from the lowest to highest mole fractions for each.

    Keyword arguments:
    interpolate_pmf_results -- dictionary output from the interpolate_pmf function
    """

    bins = []

    # Creates numpy array of all compositions, needed to calculate min/max of each gas's mole frac.
    comps_array = np.array([[row['CO2'], row['C2H6'], row['CH4'], row['N2']] for row in interpolate_pmf_results])
    bin_range = np.column_stack((np.linspace(min(comps_array[:,0]), max(comps_array[:,0]), 12),
                                 np.linspace(min(comps_array[:,1]), max(comps_array[:,1]), 12),
                                 np.linspace(min(comps_array[:,2]), max(comps_array[:,2]), 12),
                                 np.linspace(min(comps_array[:,3]), max(comps_array[:,3]), 12)))
    bins.extend([{'CO2' : row[0], 'C2H6' : row[1], 'CH4' : row[2], 'N2' : row[3]} for row in bin_range])
    return(bins)

def bin_compositions(gases,mof_array, create_bins_results, interpolate_pmf_results):
    """Sorts pmfs into bins created by create_bins function.

    Keyword arguments:
    gases -- list of gases specified as user input
    mof_array -- list of MOFs in the array, specified as user input
    create_bins_results -- dictionary containing bins for each gas
    interpolate_pmf_results -- dictionary of pmfs associated with MOFs/gases as a result of interpolate_pmf function
    """

    binned_probability = []
    for mof_name in mof_array:
        for gas_name in gases:
            binned_data =[]
            binned_probability_temporary = []

            # Assigns pmf to bin value (dictionary) by checking whether mole frac is between the current and next bin value.
            for row in interpolate_pmf_results:
                 for i in range(1, len(create_bins_results)):
                    if row[gas_name] >= create_bins_results[i - 1][gas_name] and row[gas_name] < create_bins_results[i][gas_name] and row['mof'] == mof_name:
                        binned_data.append({'probability': row['PMF 1bar'], 'bin': create_bins_results[i - 1][gas_name]})

            # Loops through all of the bins and averages the pmfs into their assgned bins.
            for b in create_bins_results:
                average = []
                for line in binned_data:
                     if b[gas_name] == line['bin']:
                        average.append(line['probability'])
                if average == []:
                    binned_probability_temporary.append({'bin' : line['bin'], 'average probability' : 0})
                else:
                    binned_probability_temporary.append({'bin' : line['bin'], 'average probability' : np.mean(average)})

            # Creates list of binned probabilities in order to loop through and normalize, sum must be 1.
            temporary_pmf_list = [row['average probability'] for row in binned_probability_temporary]
            normalized_temporary_pmf = [number / sum(temporary_pmf_list) for number in temporary_pmf_list]
            binned_probability.extend([{'mof' : mof_name, 'gas' : gas_name, 'bin' : binned_probability_temporary[i]['bin'],
                                        'average probability' : normalized_temporary_pmf[i]} for i in range(0, len(normalized_temporary_pmf))])
    return(binned_probability)

def plot_binned_pmf_array(gas_names,mof_names, bin_compositions_results, create_bins_results):
    """Calculates compound pmfs for MOF array and plots vs mole fraction for each gas.

    Keyword arguments:
    gas_names -- list of gases specified by user
    mof_names -- list of MOFs in array, specified by user
    bin_compositions_results -- dictionary result from bin_compositions function
    create_bins_results -- dictionary result from create_bins
    """

    for gas_name in gas_names:
        compound_pmfs = []
        for mof in mof_names:
            if compound_pmfs == []:
                compound_pmfs = np.array([point['average probability'] for point in bin_compositions_results if point['mof'] == mof and point['gas'] == gas_name])
            else:
                compound_pmfs *= np.array([point['average probability'] for point in bin_compositions_results if point['mof'] == mof and point['gas'] == gas_name])
        normalized_compound_pmfs = [number / sum(compound_pmfs) for number in compound_pmfs]
        plot_PMF = plt.figure()
        plt.plot([b[gas_name] for b in create_bins_results], [point for point in normalized_compound_pmfs], 'bo')
        plt.savefig("plot_PMF_%s_%s.png" % (str(gas_name) , "_".join(mof_names)))
        plt.close(plot_PMF)
