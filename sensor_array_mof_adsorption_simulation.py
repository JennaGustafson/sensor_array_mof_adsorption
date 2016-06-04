import os
import subprocess
import shutil

from jobserver_utils import generate_unique_per_process_filename

def write_raspa_file(filename, mof, co2_mf, ch4_mf, n2_mf, c2h6_mf):
    shutil.copyfile("simulation_template.file", filename)
    subprocess.run(['sed', '-i', "s/^FrameworkName .*$/FrameworkName %s/" % mof, filename], check=True)
    subprocess.run(['sed', '-i', "20 s/[0-9]*\.[0-9]*/%s/" % co2_mf,  filename], check=True)
    subprocess.run(['sed', '-i', "31 s/[0-9]*\.[0-9]*/%s/" % ch4_mf,  filename], check=True)
    subprocess.run(['sed', '-i', "42 s/[0-9]*\.[0-9]*/%s/" % n2_mf,   filename], check=True)
    subprocess.run(['sed', '-i', "53 s/[0-9]*\.[0-9]*/%s/" % c2h6_mf, filename], check=True)

def parse_output(output_file):
    mass = float(subprocess.check_output(['./calculate_mass.sh', output_file]))
    return mass

def run(mof, co2_mf, ch4_mf, n2_mf, c2h6_mf, output_dir='output'):
    # create unique working directory for this simulation
    working_dir = os.path.join(output_dir, generate_unique_per_process_filename())
    os.makedirs(working_dir, exist_ok=True)

    # run simulation
    write_raspa_file(os.path.join(working_dir, "simulation.input"), mof, co2_mf, ch4_mf, n2_mf, c2h6_mf)
    subprocess.run(['simulate', 'simulation.input'], check=True, cwd=working_dir)

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
