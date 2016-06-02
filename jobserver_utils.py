import os
from datetime import datetime

def generate_unique_run_name():
    return datetime.now().strftime("%Y_%m_%d__%H_%M_%S")

def generate_unique_per_process_filename():
    return "%s_%s" % (os.uname()[1], os.getpid())
