"""
AHRD wrapper
"""

# import libraries
import argparse
from argparse import RawTextHelpFormatter
import os
import re
import sys
import subprocess
import logging
import yaml
import glob
import pkg_resources
from shutil import copy
import pandas as pd

from eifunannot import __version__, __author__, __email__

# check python version
try:
    assert sys.version_info >= (3, 7)
except AssertionError:
    logging.error(
        f"Python >=3.7 is required.\n"
        f"Current is Python {sys.version}.\nExiting..."
    )
    sys.exit()

# check python version
if sys.version_info[0] < 3:
    raise Exception("Please source Python 3, sourcing 'source snakemake-5.4.0' will do")

# get script name
script =  os.path.basename(sys.argv[0])

# get cwd
cwd = os.getcwd()

# get pkg_resources
# get config location
# eifunannot_folder = pkg_resources.resource_filename("eifunannot", "eifunannot")
config_folder = pkg_resources.resource_filename("eifunannot", "config")
# need to split the file and make uniq list names
configure_files = [os.path.basename(fname) for fname in glob.iglob(os.path.join(config_folder, "**", "*yaml"), recursive=True)]
cluster_config = "".join([fname for fname in glob.iglob(os.path.join(config_folder, "**", "*cluster_config*json"), recursive=True)])
run_config = [fname for fname in glob.iglob(os.path.join(config_folder, "**", "*run_config*yaml"), recursive=True)]
ahrd_config = [fname for fname in glob.iglob(os.path.join(config_folder, "**", "*ahrd*generic*yaml"), recursive=True)]
# snakemake path
# snakemake_file = [fname for fname in glob.iglob(os.path.join(eifunannot_folder, "**", "eifunannot.smk"), recursive=True)]
snakemake_file = pkg_resources.resource_filename(
    "eifunannot", "eifunannot.smk")


class EiFunAnnotAHRD(object):

    def __init__(self):
        parser = argparse.ArgumentParser(
            description=f'EI FunAnnot version {__version__}',
            usage='''eifunannot <command> [<args>]

The commands are:
   configure     Get the configuration files
   run           Run the pipeline
''')
        parser.add_argument('command', help='Subcommand to run')
        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print('Unrecognized command')
            parser.print_help()
            sys.exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()

    # eifunannot configure setup
    def configure(self):

        parser = argparse.ArgumentParser(description=f"EI FunAnnot version {__version__} - configure", formatter_class=RawTextHelpFormatter, 
        epilog="Example configure command:\n" + script + " configure" + "\n\nContact:" + __author__ + "(" + __email__ + ")")

        # parser.add_argument("-c", "--configure", type=str, default=None, help="Configure file to use. Choose from the following:\n{}".format(",\n".join(configure_files)))
        parser.add_argument("-o", "--output", type=str, default=cwd, help="Output directory (default: %(default)s)")
        parser.add_argument("-f", "--force", action='store_true', help="Force overwrite if configuration files exist (default: %(default)s)")
        # now that we're inside a subcommand, ignore the first
        # TWO argvs, ie the command (eifunannot) and the subcommand (configure)
        args = parser.parse_args(sys.argv[2:])
        # configure = args.configure
        output = args.output
        force = args.force
        print(f"Running eifunannot configure..")
        for file in cluster_config, run_config, ahrd_config:
            abs_file_path = "".join(file)
            file_base = os.path.basename(abs_file_path)
            if force or not os.path.exists(os.path.join(output, file_base)):
                print(f"Copying file {file_base} to {cwd}")
                copy(abs_file_path, cwd)
            elif os.path.exists(os.path.join(output, file_base)):
                print(f"Configure file already exists {os.path.join(output, file_base)}")
                print(f"Use --force to overwrite")
                sys.exit(1)

    # eifunannot run setup
    def run(self):

        parser = argparse.ArgumentParser(description=f"EI FunAnnot version {__version__} - run", formatter_class=RawTextHelpFormatter, 
        epilog="Example run command:\n" + script + " run --config [plants.yaml]" + "\n\nContact:" + __author__ + "(" + __email__ + ")")
        parser.add_argument("--config", required=True, nargs='?',
                            help="Provide config.yaml file to run the pipeline (default: %(default)s)")
        parser.add_argument("--ahrd_config", required=True, nargs='?',
                            help="Provide ahrd_config.yaml file to run the pipeline (default: %(default)s)")
        parser.add_argument("--hpc_config", default=cluster_config,
                            nargs='?', help="Provide hpc_config.json file (default: %(default)s)")
        parser.add_argument("--jobs", type=int, default=1, nargs='?',
                            help="Maximum number of jobs to execute at any one time (default: %(default)s)")
        # boolen action='store_true' help from here: https://stackoverflow.com/a/36031646
        parser.add_argument("-np", "--dry_run",
                            action='store_true', help="Dry run (default: %(default)s)")
        # args = parser.parse_args()
        # only consider arguments after first two
        args = parser.parse_args(sys.argv[2:])
        # args = parser.parse_args()
        config = args.config
        ahrd_config = args.ahrd_config
        hpc_config = args.hpc_config
        jobs = args.jobs
        dry_run = args.dry_run

        output = None
        with open(config, 'r') as ymlfile:
            cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)
            output = os.path.abspath(cfg['output'])

        # run AHRD pipeline
        print(f"Running eifunannot run..")
        EiFunAnnotAHRD.run_ahrd(
            output, config, ahrd_config, hpc_config, dry_run, jobs)
    
    @staticmethod
    def run_ahrd(output, config, ahrd_config, hpc_config, dry_run, jobs):
        # print(output, config, hpc_config, dry_run, jobs)
        # print(type(output), type(config), type(hpc_config), type(dry_run), type(jobs))
        cmd = None
        if dry_run:
            print("Enabling dry run..")
            # print(snakemake_file)
            # sys.exit(1)
            cmd = "snakemake --snakefile " + snakemake_file \
                + " --configfile " + config \
                + " --config ahrd_config=" + ahrd_config \
                + " --cluster-config " + hpc_config \
                + " -np"
        else:
            cmd = "snakemake --latency-wait 120 --snakefile " + snakemake_file \
                + " --configfile " + config \
                + " --config ahrd_config=" + ahrd_config \
                + " --cluster-config " + hpc_config \
                + " --printshellcmds" \
                + " --jobs " + str(jobs) \
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
            unknown = data['Human-Readable-Description'].str.contains(r"^Unknown protein$")
            try:
                hits = unknown.value_counts()[False]
            except KeyError as err:
                # print(f"There are no hits that has function: {err}")
                hits = 0
            try:
                no_hits = unknown.value_counts()[True]
            except KeyError as err:
                # print(f"There are no hits that are unknown: {err}")
                no_hits = 0

            # print error if the total_data does not match hits + no_hits
            assert total_data == hits + \
                no_hits, f"ERROR: Total output {total_data} and hits {hits + no_hits} do not match"
            # get repeat association - based on two columns, 'Human-Readable-Description' and 'Interpro-ID (Description)'
            # https://stackoverflow.com/a/45709524
            repeat_associated = data['Human-Readable-Description'].str.contains(
                r"transpos|helicas") | data['Interpro-ID (Description)'].str.contains(r"transpos|helicas")
            try:
                repeat_associated_hits = repeat_associated.value_counts()[True]
            except KeyError as err:
                # print(f"There are no hits that are repeat associated: {err}")
                repeat_associated_hits = 0
            # get unknown with interproscan ids
            unknown_with_ipr = data['Human-Readable-Description'].str.contains(
                r"^Unknown protein$") & data['Interpro-ID (Description)'].str.contains(r"IPR")
            try:
                unknown_with_ipr_hits = unknown_with_ipr.value_counts()[True]
            except KeyError as err:
                # print(f"There are no unknown hits with IPR: {err}")
                unknown_with_ipr_hits = 0

            #=====#
            # print summary
            print("Overall summary:")
            print(f"Of total {total_data} input protein models:")
            print(f"{hits} ({round(hits / total_data, 4) * 100} %) - have functional annotation (of which {repeat_associated_hits} are repeat associated (transposon|transposase|helicase))")
            print(
                f"{no_hits} ({round(no_hits / total_data, 4) * 100} %) - are unknown proteins (of which {unknown_with_ipr_hits} have an interproscan id)")
            print(f"\nOutput file:\n{output_file}\n")
            #=====#

if __name__ == "__main__":
    value = EiFunAnnotAHRD()
    print(f"The return is:{value}")
