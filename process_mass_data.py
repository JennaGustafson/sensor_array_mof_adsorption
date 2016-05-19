#!/usr/bin/env python
import numpy as np
import scipy.stats as ss
from matplotlib import pyplot as plt
from math import isnan
from scipy.spatial import Delaunay
import scipy.interpolate as si
import csv
from random import random
#
def mof_density(filename):                        #function imports a tab delimited csv file as dictionary
    with open(filename,newline='') as csvfile:
        density = csv.DictReader(csvfile, delimiter="\t")
        return list(density)[0]                   #returns the first row of the file
def read_output_data(filename):                   #function imports csv file as a dictionary
    with open(filename,newline='') as csvfile:
        output_data = csv.DictReader(csvfile, delimiter="\t")
        return list(output_data)                  #returns a list contatining a list of dictionaries
def mof_names(filename):                          #reads in csvfile as a list
    with open(filename,newline='') as csvf:
        names = csv.reader(csvf,delimiter='\n')
        return list(names)
#
mof_densities_import = mof_density('MOF_Density.csv')               #employ function to import the densities of each mof, stored in a csvfile
all_results_import = read_output_data('comp_mass_output_tmp.csv')   #employ function to import simulation output data into dict format
mofs_import = [row[0] for row in mof_names('mofs.csv')]             #employ function to import mof names as a list
mof_experimental_mass = mof_density('MOF_ExperimentalMass.csv')     #employ function to import mof densities as a dictionary
#
num_mixtures = 10 #specify how many *additional* mixtures one wants in the analysis (to be interpolated)
r = .05           #the range for which the pmf is calculated, (cdf mass+r minus cdf mass-r)
stdev = 0.1       #standard deviation of the  experimental masses for normal distribution curves
#
def interpolate_pmf(mofs_list,all_results,experimental_mass,mof_densities):      #function which returns the pmfs of each mass for each mixture/mof combination
    pmf_results = []
    for mof in mofs_list:
        masses = [float(mof_densities[row['MOF']])*float(row['Mass 1bar']) for row in all_results if row['MOF'] == mof]   #miltuply mass by mof density to get mg/(cm3 of mof structure)
        comps = [[float(row['CH4']),float(row['CO2']),float(row['C2H6'])] for row in all_results if row['MOF'] == mof]    #isolate the mixture composition values
        d = Delaunay(comps)                                                      #perform delaunay tessellation on the original compositions
        interp_dat = si.LinearNDInterpolator(d,masses)                           #interpolate using the tesselation and the corresponding mass values for each mixture
        while (len(comps) < 78+num_mixtures):                                    #check that number of mixtures added does not exceed num_mixtures
            random_gas = ([0.5*round(random(),3),0.5*round(random(),3),0.2*round(random(),3)])    #generate random gas mixture, from 0-0.5 (CH4),0-0.5 (CO2), 0-0.2 (C2H6)
            predicted_mass = interp_dat(random_gas)                              #calculate the inperpolated mass at the new gas mixture
            if sum(random_gas) <= 1 and not isnan(predicted_mass):               #check that the total mole fraction does not exceed 1 and that the new mass is a real number
                comps.append(random_gas)                                         #if these are trye then append the new mixture to the original set
                masses.extend(predicted_mass)                                    #append the new mass to the original mass values
        probs = [(ss.norm.cdf(mass+r, float(experimental_mass[mof]),             #calculate probability mass functions for list of all masses
                              stdev*float(experimental_mass[mof])) -             #uses the *experimental* mass values for each mof
                  ss.norm.cdf(mass-r, float(experimental_mass[mof]),
                              stdev*float(experimental_mass[mof]))) for mass in masses]
        norm_probs = [(i/sum(probs)) for i in probs]                             #normalize the pmfs, so that they add to one
        comps_mass_prob = np.column_stack((comps,masses,norm_probs))             #construct array of composition values, masses and the pmfs
        pmf_results.extend([{'mof': mof, 'CH4': row[0], 'CO2': row[1], 'C2H6': row[2], 'N2' : 1-(row[0]+row[1]+row[2]),   #create dictionary with all info
                         'Mass 1bar' : row[3], 'PMF 1bar' : row[4]} for row in comps_mass_prob])
    return(pmf_results)                                                          #returns complete dictionary
#
def create_bins(interpolate_pmf_results):                                        #function creates arrays of bins for each gas
    bins = []
    comps_array = np.array([[row['CO2'],row['C2H6'],row['CH4'],row['N2']] for row in interpolate_pmf_results]) #separate compositions into an iterable array
    bin_range = np.column_stack((np.linspace(min(comps_array[:,0]),max(comps_array[:,0]),12),                  #uses linspace to create evenly spaced bins
                                 np.linspace(min(comps_array[:,1]),max(comps_array[:,1]),12),                  #btwn min and max of the gas's mole fracs
                                 np.linspace(min(comps_array[:,2]),max(comps_array[:,2]),12),
                                 np.linspace(min(comps_array[:,3]),max(comps_array[:,3]),12)))
    bins.extend([{'CO2' : row[0], 'C2H6' : row[1], 'CH4' : row[2], 'N2' : row[3]} for row in bin_range])       #new dictionary with the list of bins
    return(bins)
#
def bin_compositions(gases,mof_array,create_bins_results,interpolate_pmf_results):   #verb bin, sorts all pmfs into bins associated with their corresponding comps
    binned_probability = []                                                          #overall list with all binned pmfs
    for mof_name in mof_array:
        for gas_name in gases:                                                       #iterates through every mof/gas combo specified in input argumants
            binned_data =[]                                                          #dictionary with all points and their associated bins
            binned_probability_temporary = []                                        #dictionary which is overwritten for each mof/gas combo
            for row in interpolate_pmf_results:
                 for i in range(1,len(create_bins_results)):                         #sorts pmfs into bins by checking whether or not they are between the current bin value and bin value which is one ahead
                    if row[gas_name]>=create_bins_results[i-1][gas_name] and row[gas_name]<create_bins_results[i][gas_name] and row['mof'] == mof_name:
                        binned_data.append({'probability': row['PMF 1bar'], 'bin': create_bins_results[i-1][gas_name]})
            for b in create_bins_results:                                            #for each of the bin values in the bins
                average = []
                for line in binned_data:                                             #for each data point which has been binned
                     if b[gas_name] == line['bin']:                                  #check if the data point's bin matches the current bin value
                        average.append(line['probability'])                          #if they match, append the list with the pmf value
                if average == []:                                                    #if there are no points which fall into the current bin
                    binned_probability_temporary.append({'bin' : line['bin'], 'average probability' : 0})     #set the pmf bin value to be zero
                else:
                    binned_probability_temporary.append({'bin' : line['bin'], 'average probability' : np.mean(average)})  #otherwise take the average of all pmfs in the current bin
            temporary_pmf_list = [row['average probability'] for row in binned_probability_temporary]         #new list of all pmfs after binning
            normalized_temporary_pmf = [number/sum(temporary_pmf_list) for number in temporary_pmf_list]      #normalize so that total pmf over all bins is 1
            binned_probability.extend([{'mof' : mof_name, 'gas' : gas_name, 'bin' : binned_probability_temporary[i]['bin'],    #new dictionary includes all info
                    'average probability' : normalized_temporary_pmf[i]} for i in range(0,len(normalized_temporary_pmf))])     #for each bin with averaged/normalized pmf
    return(binned_probability)
#
def plot_binned_pmf_array(gas_names,mof_names,bin_compositions_results,create_bins_results):    #function compounds pmfs for mof array and plots results each array/gas
    for gas_name in gas_names:
        compound_pmfs = []
        for mof in mof_names:
            if compound_pmfs == []:                                                   #multiplies pmfs for mofs included in the array
                compound_pmfs = np.array([point['average probability'] for point in bin_compositions_results if point['mof'] == mof and point['gas'] == gas_name])
            else:
                compound_pmfs *= np.array([point['average probability'] for point in bin_compositions_results if point['mof'] == mof and point['gas'] == gas_name])
        normalized_compound_pmfs = [number/sum(compound_pmfs) for number in compound_pmfs]      #normalizes the compound probabilities
        plot_PMF = plt.figure()
        plt.plot([b[gas_name] for b in create_bins_results],[point for point in normalized_compound_pmfs],'bo')  #plot the compound pmfs vs comp bins for each mof/gas
        plt.savefig("plot_PMF_%s_%s.png" % (str(gas_name) , "_".join(mof_names)))                                #save the plot with specific figure name
        plt.close(plot_PMF)

gases = ['N2','CH4','CO2','C2H6']
mof_array = ['IRMOF-1','HKUST-1','NU-125']
interpolate_pmf_results = interpolate_pmf(mofs_import,all_results_import,mof_experimental_mass,mof_densities_import)
create_bins_results = create_bins(interpolate_pmf_results)
bin_compositions_results = bin_compositions(gases,mof_array,create_bins_results,interpolate_pmf_results)
plot_binned_pmf_array(gases,mof_array,bin_compositions_results,create_bins_results)
