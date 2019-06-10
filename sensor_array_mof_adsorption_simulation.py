import os
import subprocess
import shutil

import yaml

from jobserver_utils import generate_unique_per_process_filename

def yaml_loader(filepath):
    with open(filepath, 'r') as yaml_file:
        data = yaml.load(yaml_file)
    return(data)

def write_raspa_file(filename, mof, pressure, gases, composition, config_file):
    config_data = yaml_loader(config_file)
    gas_names_def = config_data['Forcefield_Gas_Names']
    f = open(filename,'w',newline='')

    simulation_file_header = """
SimulationType                MonteCarlo
NumberOfCycles                2000
NumberOfInitializationCycles  1000
PrintEvery                    200

ChargeMethod                  Ewald
CutOff                        12.0
Forcefield                    JennaUFF2
EwaldPrecision                1e-6

Framework 0
FrameworkName %s
UnitCells 1 1 1
HeliumVoidFraction 0.81
UseChargesFromCIFFile yes
ExternalTemperature 298.0
ExternalPressure %s
""" % (mof, pressure)

    f.write(simulation_file_header)

    component_number = 0
    for gas in gases:
        gas_name = gas_names_def[gas]
        mole_fraction = composition[gas]
        simulation_file_gas = """
    Component %s MoleculeName              %s
                 MoleculeDefinition         TraPPE-Zhang
                 MolFraction                %s
                 TranslationProbability     0.5
                 RegrowProbability          0.5
                 IdentityChangeProbability  1.0
                   NumberOfIdentityChanges  2
                   IdentityChangesList      0 1
                 SwapProbability            1.0
                 CreateNumberOfMolecules    0

                 """ % (component_number, gas_name, mole_fraction)

        f.write(simulation_file_gas)
        component_number += 1

    f.close()

def parse_output(output_file):
    mass = float(subprocess.check_output(['./calculate_mass.sh', output_file]))
    return mass

def run(run_id, mof, pressure, gases, composition, config_file, output_dir='output'):
    # create unique working directory for this simulation
    working_dir = os.path.join(output_dir, generate_unique_per_process_filename())
    os.makedirs(working_dir, exist_ok=True)

    # run simulation
    write_raspa_file(os.path.join(working_dir, "simulation.input"), mof, pressure, gases,
                        composition, config_file
                    )
    subprocess.run(['simulate', 'simulation.input'], check=True, cwd=working_dir)

    # parse data from simulation
    data_filename = os.path.join(working_dir, 'Output', 'System_0', '*.data')
    mass = parse_output(data_filename)

    # archive data and configuration; delete working_dir
    run_descriptor = "%s" % (run_id)
    archive_dir = os.path.join(output_dir, 'archive', run_descriptor)
    os.makedirs(archive_dir, exist_ok=True)
    for f in ["Output", "simulation.input"]:
        shutil.move(os.path.join(working_dir, f), archive_dir)

    shutil.rmtree(os.path.join(working_dir))

    return (mass)
