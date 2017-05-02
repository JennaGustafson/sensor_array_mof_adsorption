#!/usr/bin/env python3
import os
import sys
import csv

from jobserver_utils import generate_unique_run_name
import sjs
from datetime import datetime
from mof_array.sim.sensor_array_mof_adsorption import (read_composition_configuration,
                                                      read_mof_configuration,
                                                      run_composition_simulation)

mofs_filepath = sys.argv[1]
gas_comps_filepath = sys.argv[2]
pressure = sys.argv[3]

mofs = read_mof_configuration(mofs_filepath)
compositions = read_composition_configuration(gas_comps_filepath)
gases = list(compositions[0].keys())

run_name = generate_unique_run_name()
output_dir = 'output_%s' % run_name
os.makedirs(output_dir)

sjs.load(os.path.join("settings","sjs.yaml"))
job_queue = sjs.get_job_queue()

# setup CSV file and write header
f = open(os.path.join(output_dir, 'comp_mass_output.csv'),'w',newline='')
# write header
header = ['Run ID','MOF','Mass']
for gas in gases:
    header.append(gas)
writer = csv.writer(f, delimiter='\t')
writer.writerow(header)

if job_queue is not None:
    print("Queueing jobs onto queue: %s" % job_queue)

    run_id_number = 0
    for mof in mofs:
        for composition in compositions:
            run_id = "%s" % run_id_number
            job_queue.enqueue(run_composition_simulation, run_id, mof, pressure, gases, composition, csv_writer=None, output_dir=output_dir)
            run_id_number += 1

else:
    print("No job queue is setup. Running in serial mode here rather than on the cluster")

    run_id_number = 0
    for mof in mofs:
        for composition in compositions:
            run_id = "%s_%s" % (run_id_number, datetime.now().strftime("%Y_%m_%d__%H_%M_%S"))
            run_composition_simulation(run_id, mof, pressure, gases, composition, csv_writer=writer, output_dir=output_dir)
            run_id_number +=1

f.close()
