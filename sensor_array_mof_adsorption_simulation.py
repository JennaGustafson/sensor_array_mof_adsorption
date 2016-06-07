import os
import subprocess
import shutil

import yaml

def yaml_loader(filepath):
    with open(filepath, 'r') as yaml_file:
        data = yaml.load(yaml_file)
    return(data)

def write_raspa_file(mof, pressure, gases, composition, config_file):
    config_data = yaml_loader(config_file)
    gas_names_def = config_data['Forcefield_Gas_Names']
    f = open('simulate.input','w',newline='')

    simulation_file_header = """
SimulationType                MonteCarlo
NumberOfCycles                1000
NumberOfInitializationCycles  1000
PrintEvery                    100

ChargeMethod                  Ewald
CutOff                        12.0
Forcefield                    GenericMOFs
EwaldPrecision                1e-6

Framework 0
FrameworkName %s
UnitCells 1 1 1
HeliumVoidFraction 0.81
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
                 MoleculeDefinition         TraPPE
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

def run(run_id, mof, pressure, gases, composition, config_file):
    write_raspa_file(mof, pressure, gases, composition, config_file)
    subprocess.run(['simulate', 'simulate.input'], check=True)

    mass = parse_output('./Output/System_0/*00.data')

    # archive data and configuration
    run_descriptor = "%s" % (run_id)
    archive_path = os.path.join('archive', run_descriptor)
    os.makedirs(archive_path, exist_ok=True)
    shutil.move('simulate.input', archive_path)
    shutil.move('Output', archive_path)

    # delete extra directories
    shutil.rmtree('Movies')
    shutil.rmtree('VTK')
    shutil.rmtree('Restart')

    return (mass)
