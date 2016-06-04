import os
import subprocess
import shutil

from jobserver_utils import generate_unique_per_process_filename

def write_raspa_file(filename, mof, pressure, gases, composition):
    f = open(filename,'w',newline='')

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
""" $(mof, pressure)

    f.write(simulation_file_header)

    component_number = 0
    for gas in gases:
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

                 """ $(component_number, gas, mole_fraction)

        f.write(simulation_file_gas)
        component_number += 1

    f.close()

def parse_output(output_file):
    mass = float(subprocess.check_output(['./calculate_mass.sh', output_file]))
    return mass

def run(run_descriptor, mof, pressure, gases, composition, output_dir='output'):
    # create unique working directory for this simulation
    working_dir = os.path.join(output_dir, generate_unique_per_process_filename())
    os.makedirs(working_dir, exist_ok=True, cwd=working_dir)

    # run simulation
    write_raspa_file(os.path.join(working_dir, "simulation.input"), mof, pressure, gases, composition)
    subprocess.run(['simulate', 'simulate.input'], check=True)

    # parse data from simulation
    data_filename = os.path.join(working_dir, 'Output', 'System_0', '*00.data')
    mass = parse_output(data_filename)

    # archive data and configuration; delete working_dir
    run_descriptor = "%s_%s_%s_%s_%s" % (mof, co2_mf, ch4_mf, n2_mf, c2h6_mf)
    archive_dir = os.path.join(output_dir, 'archive', run_descriptor)
    os.makedirs(archive_dir, exist_ok=True)
    for f in ["Output", "simulation.input"]:
        shutil.move(os.path.join(working_dir, f), archive_dir)

    shutil.rmtree(os.path.join(working_dir))

    return (mass)
