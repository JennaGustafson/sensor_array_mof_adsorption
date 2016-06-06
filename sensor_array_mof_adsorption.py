
import os
import csv

import sensor_array_mof_adsorption_simulation

def generate_unique_per_process_filename(extension='csv'):
    pid = os.getpid()
    hostname = os.uname()[1]
    return "%s_%s.%s" % (pid, hostname, extension)

def read_mof_configuration(filename):
    with open(filename) as f:
        return [ line.strip() for line in f.readlines() ]

def read_composition_configuration(filename):
    with open(filename,newline='') as csvfile:
        comp_reader = csv.DictReader(csvfile, delimiter="\t")
        return list(comp_reader)

def run_composition_simulation(mof, composition, csv_writer=None):

    # if there is no csv_writer passed, we write to a file that is unique to this process
    csv_file = None
    if csv_writer is None:
        filename = generate_unique_per_process_filename()
        csv_file = open(filename,'a',newline='')
        csv_writer = csv.writer(csv_file, delimiter='\t')

    # run the simulation / output the data
    mass_p1, mass_p2 = sensor_array_mof_adsorption_simulation.run(
        mof,
        composition['CO2'], composition['CH4'], composition['N2'], composition['C2H6']
    )

    csv_writer.writerow([
        mof,
        composition['CO2'], composition['CH4'], composition['N2'], composition['C2H6'],
        mass_p1, mass_p2
    ])

    # close the file, if we opened it above
    if csv_file is not None:
        csv_file.close()
