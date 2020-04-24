#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# Tuesday, 28 January 2020, 10:55AM
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
import logging
import yaml
import pandas as pd

# check python version
try:
    assert sys.version_info >= (3, 7)
except AssertionError:
    logging.error(
        f"Python >=3.7 is required.\n"
        f"Current is Python {sys.version}.\nExiting..."
    )
    sys.exit(1)

# check python version
if sys.version_info[0] < 3:
    raise Exception("Please source Python 3, sourcing 'source snakemake-5.4.0' will do")

# get script name
script =  os.path.basename(sys.argv[0])

# snakemake path
snakemake_file = "/ei/cb/common/Scripts/eifunannot/0.2/AHRD_pipeline.Snakefile.smk"

def run_ahrd(output, config, hpc_config, dry_run, jobs):
    cmd = None
    if dry_run:
        print("Enabling dry run..")
        cmd = "snakemake --snakefile " + snakemake_file \
              + " --configfile " + config \
              + " --cluster-config " + hpc_config \
              + " -np"
    else:
        cmd = "snakemake --latency-wait 120 --snakefile " + snakemake_file \
              + " --configfile " + config \
              + " --cluster-config " + hpc_config \
              + " --printshellcmds" \
              + " --jobs " + jobs \
              + " --cluster \" sbatch -p {cluster.partition} -c {cluster.c} --mem {cluster.mem} -J {cluster.J} -o {cluster.o}\""
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,universal_newlines=True) # for universal_newlines - https://stackoverflow.com/a/4417735
    (result, error) = p.communicate() # https://www.programcreek.com/python/example/50/subprocess.Popen
    exit_code = p.returncode
    if exit_code:
        raise subprocess.CalledProcessError(exit_code, cmd)
    print (result)

    if dry_run:
        print("Dry run completed successfully!\n")
    else:
        print("AHRD pipeline completed successfully!\n")
        output_file = f"{output}/ahrd_output.csv"
        data = pd.read_csv(output_file, sep='\t', comment='#')
        # keeping for reference
        # print(data['AHRD-Quality-Code'].value_counts(dropna=False))
        # print(data['AHRD-Quality-Code'].value_counts(normalize=True, dropna=False)*100)
        # # get quality code that are not NaN
        # hits = data['AHRD-Quality-Code'].isna().value_counts()[False]
        # # get quality code that are NaN
        # no_hits = data['AHRD-Quality-Code'].isna().value_counts()[True]

        #=====#
        # print stats
        # AHRD-Quality-Code stats
        # https://stackoverflow.com/a/50169502
        ahrd_quality_code_stats = pd.concat([data['AHRD-Quality-Code'].value_counts(dropna=False), round(data['AHRD-Quality-Code'].value_counts(
            normalize=True, dropna=False).mul(100), 2)], axis=1, keys=('counts', 'percentages'))
        print(
            f"AHRD-Quality-Code stats: [NaN = No hit, Reference here: https://github.com/groupschoof/AHRD/blob/master/README.textile#241-tab-delimited-table]\n{ahrd_quality_code_stats}\n")
        #=====#
        total_data = len(data)  # get total
        # get count with and without function
        unknown = data['Human-Readable-Description'].str.contains(
            r"^Unknown protein$")
        hits = unknown.value_counts()[False]
        no_hits = unknown.value_counts()[True]
        # print error if the total_data does not match hits + no_hits
        assert total_data == hits + \
            no_hits, f"ERROR: Total output {total_data} and hits {hits + no_hits} do not match"
        # get repeat association - based on two columns, 'Human-Readable-Description' and 'Interpro-ID (Description)'
        # https://stackoverflow.com/a/45709524
        repeat_associated = data['Human-Readable-Description'].str.contains(
            r"transpos|helicas") | data['Interpro-ID (Description)'].str.contains(r"transpos|helicas")
        repeat_associated_hits = repeat_associated.value_counts()[True]
        # get unknown with interproscan ids
        unknown_with_ipr = data['Human-Readable-Description'].str.contains(
            r"^Unknown protein$") & data['Interpro-ID (Description)'].str.contains(r"IPR")
        unknown_with_ipr_hits = unknown_with_ipr.value_counts()[True]

        #=====#
        # print summary
        print("Overall summary:")
        print(f"Of total {total_data} input protein models:")
        print(f"{hits} ({round(hits / total_data, 4) * 100} %) - have functional annotation (of which {repeat_associated_hits} are repeat associated (transposon|transposase|helicase))")
        print(
            f"{no_hits} ({round(no_hits / total_data, 4) * 100} %) - are unknown proteins (of which {unknown_with_ipr_hits} have an interproscan id)")
        print(f"\nOutput file:\n{output_file}\n")
        #=====#


def main():
    parser = argparse.ArgumentParser(description="AHRD wrapper", formatter_class=RawTextHelpFormatter,
                                     epilog="Example command:\n"
                                     + script + " --config [config.yaml]"
                                     + " --hpc_config [hpc_config.json]"
                                     + " --jobs 100" + "\n\nContact:" + __author__ + "(" + __email__ + ")")
    parser.add_argument("--config", required=True, nargs='?',
                        help="Provide config.yaml file to run the pipeline")
    parser.add_argument("--hpc_config", required=True,
                        nargs='?', help="Provide hpc_config.json file")
    parser.add_argument("--jobs", required=True, nargs='?',
                        help="Maximum number of jobs to execute at any one time")
    # boolen action='store_true' help from here: https://stackoverflow.com/a/36031646
    parser.add_argument("-np", "--dry_run",
                        action='store_true', help="Dry run")
    args = parser.parse_args()
    config = args.config
    hpc_config = args.hpc_config
    jobs = args.jobs
    dry_run = args.dry_run

    output = None
    with open(config, 'r') as ymlfile:
        cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)
        output = os.path.abspath(cfg['output'])

    # run AHRD pipeline
    run_ahrd(output, config, hpc_config, dry_run, jobs)

if __name__ == "__main__":
    main()
