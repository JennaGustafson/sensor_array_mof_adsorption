#!/usr/bin/env python

#read in csv file with gas compositions
#read in csv file with mofs
import csv
import subprocess

def read_mof_configuration(filename):
    with open(filename) as f:
        return [ line.strip() for line in f.readlines() ]

def read_composition_configuration(filename):
    with open(filename,newline='') as csvfile:
        comp_reader = csv.DictReader(csvfile, delimiter="\t")
        return list(comp_reader)

compositions = read_composition_configuration('comps.csv')
mofs = read_mof_configuration('mofs.csv')

print(compositions)
print(mofs)

with open('comp_mass_output.csv','w',newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['MOF','CO2','CH4','N2','C2H6','Mass 1bar','Mass 10bar'])

for mof in mofs:
    for composition in compositions:
        subprocess.run(["./run_adsorption.sh",composition['CO2'],composition['CH4'],composition['N2'],composition['C2H6'],mof])
