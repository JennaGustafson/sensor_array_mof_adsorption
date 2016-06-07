#!/usr/bin/env python

#read in csv file with gas compositions
#read in csv file with mofs
import sys
import csv
import subprocess

import sensor_array_mof_adsorption_simulation

def read_mof_configuration(filename):
    with open(filename) as f:
        return [ line.strip() for line in f.readlines() ]

def read_composition_configuration(filename):
    with open(filename,newline='') as csvfile:
        comp_reader = csv.DictReader(csvfile, delimiter="\t")
        return list(comp_reader)


mofs_filepath = sys.argv[1]
gas_comps_filepath = sys.argv[2]
gases_filepath = sys.argv[3]
pressure = sys.argv[4]

compositions = read_composition_configuration(gas_comps_filepath)
mofs = read_mof_configuration(mofs_filepath)
gases = read_mof_configuration(gases_filepath)

print(compositions)
print(mofs)

f = open('comp_mass_output.csv','w',newline='')

# write header
header = ['Run ID','MOF','Mass']
for gas in gases:
    header.append(gas)

writer = csv.writer(f, delimiter='\t')
writer.writerow(header)

run_id = 67161308

for mof in mofs:
    run__id = run_id + 111
    for composition in compositions:
        mass = sensor_array_mof_adsorption_simulation.run(
            run_id, mof, pressure, gases, composition, 'write_comps_config.yaml'
        )

        writer.writerow([
            run_id, mof, mass, "	".join( str(composition[gas]) for gas in gases)
        ])
        run_id += 1

f.close()
