#!/usr/bin/env python

import numpy as np
import scipy.stats as ss
from matplotlib import pyplot as plt
from math import isnan
from scipy.spatial import Delaunay
import scipy.interpolate as si
import csv
from random import random
plt.close('all')
#
def mof_density(filename):
    with open(filename,newline='') as csvfile:
        density = csv.DictReader(csvfile, delimiter="\t")
        return list(density)[0]
def read_output_data(filename):
    with open(filename,newline='') as csvfile:
        output_data = csv.DictReader(csvfile, delimiter="\t")
        return list(output_data)
def mof_names(filename):
    with open(filename,newline='') as csvf:
        names = csv.reader(csvf,delimiter='\n')
        return list(names)
#
mof_densities = mof_density('MOF_Density.csv')
all_results = read_output_data('comp_mass_output_tmp.csv')
mofs = [row[0] for row in mof_names('mofs.csv')]
mof_experimental_mass = mof_density('MOF_ExperimentalMass.csv')
#
results = []
num_mixtures = 10
r = .05
stdev = 0.1
#
def interpolate_pmf(mofs,results,experimental_mass,mof_densities):
    for mof in mofs:
        print(mof)
        masses = [float(mof_densities[row['MOF']])*float(row['Mass 1bar']) for row in results if row['MOF'] == mof]
        comps = [[float(row['CH4']),float(row['CO2']),float(row['C2H6'])] for row in results if row['MOF'] == mof]
        d = Delaunay(comps)
        interp_dat = si.LinearNDInterpolator(d,masses)
        while (len(comps) < 78+num_mixtures):
            random_gas = ([0.5*round(random(),3),0.5*round(random(),3),0.2*round(random(),3)])
            predicted_mass = interp_dat(random_gas)
            if sum(random_gas) <= 1 and not isnan(predicted_mass):
                comps.append(random_gas)
                masses.extend(predicted_mass)
        probs = [(ss.norm.cdf(mass+r, float(experimental_mass[mof]),
                              stdev*float(experimental_mass[mof])) -
                  ss.norm.cdf(mass-r, float(experimental_mass[mof]),
                              stdev*float(experimental_mass[mof]))) for mass in masses]
        norm_probs = [(i/sum(probs)) for i in probs]
        comps_mass_prob = np.column_stack((comps,masses,norm_probs))
        results.extend([{'mof': mof, 'CH4': row[0], 'CO2': row[1], 'C2H6': row[2], 'N2' : 1-(row[0]+row[1]+row[2]),
                         'Mass 1bar' : row[3], 'PMF 1bar' : row[4]} for row in comps_mass_prob])
    return(results)
#
def create_bins(interpolate_pmf_results):
    bins = []
    comps_array = np.array([[i['CO2'],i['C2H6'],i['CH4'],i['N2']] for i in interpolate_pmf_results])
    bin_range = np.column_stack((np.linspace(min(comps_array[:,0]),max(comps_array[:,0]),12),
                                 np.linspace(min(comps_array[:,1]),max(comps_array[:,1]),12),
                                 np.linspace(min(comps_array[:,2]),max(comps_array[:,2]),12),
                                 np.linspace(min(comps_array[:,3]),max(comps_array[:,3]),12)))
    bins.extend([{'CO2' : row[0], 'C2H6' : row[1], 'CH4' : row[2], 'N2' : row[3]} for row in bin_range])
    return(bins)
#
def bin_compositions(gases,mof_array,create_bins_results,interpolate_pmf_results):
    binned_data =[]
    binned_probability = []
    for mof_name in mof_array:
        for gas_name in gases:
            for row in interpolate_pmf_results:
                 for i in range(1,len(create_bins_results)):
                    if row[gas_name]>=create_bins_results[i-1][gas_name] and row[gas_name]<create_bins_results[i][gas_name]:
                        binned_data.append({'probability': row['PMF 1bar'], 'bin': create_bins_results[i-1][gas_name]})
            for b in create_bins_results:
                average = []
                for line in binned_data:
                     if b[gas_name] == line['bin']:
                        average.append(line['probability'])
                if average == []:
                    binned_probability.append({'mof' : mof, 'gas' : gas_name, 'bin' : line['bin'], 'average probability' : 0})
                else:
                    binned_probability.append({'mof' : mof, 'gas' : gas_name, 'bin' : line['bin'], 'average probability' : np.mean(average)})
    return(binned_probability)
#
def plot_binned_pmf(gas_name,mof_name):
    bins_new = bin_compositions(gas_name,mof_name)
    plot_PMF = plt.figure()
    plt.plot([b[gas_name] for b in bins],[point['average probability'] for point in bins_new],'ro')
    plt.savefig("plot_PMF_%s_%s.png" % (str(gas_name) , str(mof_name)))
    plt.close(plot_PMF)
#
def plot_binned_pmf_array(gas_names,mof_names,bin_compositions_results,create_bins_results):
    for gas_name in gas_names:
        compound_pmfs = []
        for mof in mof_names:
            if compound_pmfs == []:
                compound_pmfs = np.array([point['average probability'] for point in bin_compositions_results if point['mof'] == mof and point['gas'] == gas_name])
            else:
                compound_pmfs *= np.array([point['average probability'] for point in bin_compositions_results if point['mof'] == mof and point['gas'] == gas_name])
            print(compound_pmfs)
        print(len(compound_pmfs))
        plot_PMF = plt.figure()
        plt.plot([b[gas_name] for b in create_bins_results],[point for point in compound_pmfs],'bo')
        plt.savefig("plot_PMF_%s_%s.png" % (str(gas_name) , "_".join(mof_names)))
        plt.close(plot_PMF)

gases = ['N2','CH4','CO2','C2H6']
mof_array = ['IRMOF-1','HKUST-1']
interpolate_pmf_results = interpolate_pmf(mofs,all_results,mof_experimental_mass,mof_densities)
create_bins_results = create_bins(interpolate_pmf_results)
bin_compositions_results = bin_compositions(gases,mof_array,create_bins_results,interpolate_pmf_results)
plot_binned_pmf_array(gases,mof_array,bin_compositions_results,create_bins_results)
