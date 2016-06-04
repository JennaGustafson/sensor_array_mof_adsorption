import os
import subprocess
import shutil

def write_raspa_file(filename, mof, co2_mf, ch4_mf, n2_mf, c2h6_mf):
    shutil.copyfile("simulation_template.file",filename)
    subprocess.run(['sed', '-i', "s/^FrameworkName .*$/FrameworkName %s/" % mof, filename], check=True)
    subprocess.run(['sed', '-i', "20 s/[0-9]*\.[0-9]*/%s/" % co2_mf,  filename], check=True)
    subprocess.run(['sed', '-i', "31 s/[0-9]*\.[0-9]*/%s/" % ch4_mf,  filename], check=True)
    subprocess.run(['sed', '-i', "42 s/[0-9]*\.[0-9]*/%s/" % n2_mf,   filename], check=True)
    subprocess.run(['sed', '-i', "53 s/[0-9]*\.[0-9]*/%s/" % c2h6_mf, filename], check=True)

def parse_output(output_file):
    mass = float(subprocess.check_output(['./calculate_mass.sh', output_file]))
    return mass

def run(mof, co2_mf, ch4_mf, n2_mf, c2h6_mf):
    write_raspa_file('simulate.input', mof, co2_mf, ch4_mf, n2_mf, c2h6_mf)
    subprocess.run(['simulate', 'simulate.input'], check=True)

    mass_p1 = parse_output('./Output/System_0/*00.data')

    # archive data and configuration
    run_descriptor = "%s_%s_%s_%s_%s" % (mof, co2_mf, ch4_mf, n2_mf, c2h6_mf)
    archive_path = os.path.join('archive', run_descriptor)
    os.makedirs(archive_path, exist_ok=True)
    shutil.move('simulate.input', archive_path)
    shutil.move('Output', archive_path)

    # delete extra directories
    shutil.rmtree('Movies')
    shutil.rmtree('VTK')
    shutil.rmtree('Restart')

    return (mass)
