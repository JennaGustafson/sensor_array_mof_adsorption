#!/usr/bin/env python3
import sys
import csv

from sensor_array_mof_adsorption import read_composition_configuration, read_mof_configuration
from sensor_array_mof_adsorption import run_composition_simulation
from hpc import job_queue

mofs_filepath = sys.argv[1]
gas_comps_filepath = sys.argv[2]

mofs = read_mof_configuration(mofs_filepath)
compositions = read_composition_configuration(gas_comps_filepath)

if job_queue is not None:
    print("Queueing jobs onto queue: %s" % job_queue)
    for mof in mofs:
        for composition in compositions:
            job_queue.enqueue(run_composition_simulation, mof, composition, None)

else:
    print("No job queue is setup. Running in serial mode here rather than on the cluster")

    # setup CSV file and write header
    f = open('comp_mass_output.csv','w',newline='')
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['MOF','CO2','CH4','N2','C2H6','Mass 1bar','Mass 10bar'])

    for mof in mofs:
        for composition in compositions:
            run_composition_simulation(mof, composition, writer)

    f.close()
