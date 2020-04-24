#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
05/04/2019, 11:44:00
AHRD wrapper
"""

# authorship and License information
__author__ = "Gemy George Kaithakottil"
__copyright__ = "Copyright 2018"
__license__ = "GNU General Public License v3.0"
__maintainer__ = "Gemy George Kaithakottil"
__email__ = "Gemy.Kaithakottil@gmail.com"
__status__ = "Production"
__version__ = "0.1"


# import libraries
import argparse
from argparse import RawTextHelpFormatter
import os
import re
import sys
import subprocess
# from enum import Enum, unique

# check python version
if sys.version_info[0] < 3:
    raise Exception("Please source Python 3, sourcing 'source snakemake-5.4.0' will do")

# get script name
script =  os.path.basename(sys.argv[0])


def main():
    parser = argparse.ArgumentParser(description="AHRD wrapper", formatter_class=RawTextHelpFormatter,
            epilog="Example command:\n"
            + script + " --config [config.yaml]"
            + " --hpc_config [hpc_config.json]"
            + " --jobs 100" + "\n\nContact:" + __author__ + "(" + __email__ + ")")
    parser.add_argument("--config", required=True, nargs='?', help="Provide config.yaml file to run the pipeline")
    parser.add_argument("--hpc_config", required=True, nargs='?', help="Provide hpc_config.json file")
    parser.add_argument("--jobs", required=True, nargs='?', help="Maximum number of jobs to execute at any one time")
    parser.add_argument("-np","--dry_run", action='store_true', help="Dry run") # boolen action='store_true' help from here: https://stackoverflow.com/a/36031646
    args = parser.parse_args()
    config = args.config
    hpc_config = args.hpc_config
    jobs = args.jobs
    dry_run = args.dry_run

    # run AHRD pipeline
    run_ahrd(config,hpc_config,dry_run,jobs)

def run_ahrd(config,hpc_config,dry_run,jobs):
    # cmd_test = "lsl"
    cmd = None
    if dry_run:
        cmd = "snakemake --snakefile /hpc-home/kaithakg/snakemake_scripts/AHRD_pipeline/0.1/AHRD_pipeline.Snakefile.smk" \
              + " --configfile " + config \
              + " --cluster-config " + hpc_config \
              + " -np"

    else:
        # print("No command will execute")
        cmd = "snakemake --latency-wait 120 --snakefile /hpc-home/kaithakg/snakemake_scripts/AHRD_pipeline/0.1/AHRD_pipeline.Snakefile.smk" \
              + " --configfile " + config \
              + " --cluster-config " + hpc_config \
              + " --printshellcmds" \
              + " --jobs " + jobs \
              + " --cluster \" sbatch -p {cluster.partition} -c {cluster.c} --mem {cluster.mem} -J {cluster.J} -o {cluster.o}\""
    # else:
    #     # print("No command will execute")
    #     cmd = "snakemake --snakefile /hpc-home/kaithakg/snakemake_scripts/ahrd.rodents.generic.smk" \
    #           + " --configfile " + config \
    #           + " --cluster-config " + hpc_config \

    # cmd_list = ["snakemake","--snakefile","/hpc-home/kaithakg/snakemake_scripts/ahrd.rodents.generic.smk","--configfile",config,"--cluster-config",hpc_config," -np"]
    # process = subprocess.run(cmd, returncode=0) # https://docs.python.org/3/library/subprocess.html
    # process.wait()
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,universal_newlines=True) # for universal_newlines - https://stackoverflow.com/a/4417735
    # p = subprocess.Popen(cmd_list, stdout=subprocess.PIPE,universal_newlines=True)
    (result, error) = p.communicate() # https://www.programcreek.com/python/example/50/subprocess.Popen
    exit_code = p.returncode
    # print (exit_code)
    if exit_code:
        raise subprocess.CalledProcessError(exit_code, cmd)
    print (result)

    # p_stdout = p.communicate()[0]
    # # p_stderr = p.communicate()[1]
    # print (p_stdout)
    # print (p_stderr)



if __name__ == "__main__":
    main()
