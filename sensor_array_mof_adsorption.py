import os
import csv

import sensor_array_mof_adsorption_simulation
from jobserver_utils import generate_unique_per_process_filename


def read_mof_configuration(filename):
    with open(filename) as f:
        return [ line.strip() for line in f.readlines() ]

def read_composition_configuration(filename):
    with open(filename,newline='') as csvfile:
        comp_reader = csv.DictReader(csvfile, delimiter="\t")
        return list(comp_reader)

def run_composition_simulation(mof, composition, csv_writer=None, output_dir='output'):
    # if there is no csv_writer passed, we write to a file that is unique to this process
    csv_file = None
    if csv_writer is None:
        results_dir = os.path.join(output_dir,'results')
        os.makedirs(results_dir, exist_ok=True)
        filename = os.path.join(results_dir, generate_unique_per_process_filename() + ".csv")
        csv_file = open(filename,'a',newline='')
        csv_writer = csv.writer(csv_file, delimiter='\t')

    # run the simulation / output the data
    mass = sensor_array_mof_adsorption_simulation.run(
        mof,
        composition['CO2'], composition['CH4'], composition['N2'], composition['C2H6'],
        output_dir=output_dir
    )

    csv_writer.writerow([
        mof,
        composition['CO2'], composition['CH4'], composition['N2'], composition['C2H6'],
        mass
    ])

    # close the file, if we opened it above
    if csv_file is not None:
        csv_file.close()
