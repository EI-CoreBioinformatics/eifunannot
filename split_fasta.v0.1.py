#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to split fasta file into user defined chunks


14/03/2019, 12:22:10
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
import logging
import glob
from itertools import islice
# check python version
if sys.version_info[0] < 3:
    raise Exception("Please source Python 3, sourcing 'source snakemake-5.4.0' will do")

# get working directory
cwd = os.getcwd()
# get script name
script =  os.path.basename(sys.argv[0])

# get the GTF/GFF3 attributes
SEQID, SOURCE, TYPE, START, END, SCORE, STRAND, PHASE, ATTRIBUTE  = range(9)




def main():
    parser = argparse.ArgumentParser(description="Script to split fasta file into user defined chunks", formatter_class=RawTextHelpFormatter,
            epilog="Example command:\n\t"
            + script + " --file [file.fa]"
            "\n\nContact:" + __author__ + "(" + __email__ + ")")
    parser.add_argument("-f","--file", required=True, nargs='?', help="Provide input FASTA file")
    parser.add_argument("-p","--prefix", default="chunk", nargs='?', help="Prefix for output file names [Default = \"chunk\"]")
    parser.add_argument("-c","--count", default=1000, nargs='?', type=int, help="Count of fasta to be in each chunk [Default = 1000]")
    parser.add_argument("-o","--output_dir", default=cwd,  nargs='?', help="Output directory path.\n[Default:" + cwd + "]")
    parser.add_argument("-l","--lines", type=int, default=20, help="Number of lines checked to confirm fasta format [Default = 20]")
    parser.add_argument("-v","--verbose", action="store_const", dest="loglevel", const=logging.INFO, help="Verbose output, [logging.INFO] level")
    parser.add_argument("-d","--debug", action="store_const", dest="loglevel", const=logging.DEBUG,
    default=logging.WARNING, help="Debugging messages, [logging.{WARN,DEBUG}] level")
    args = parser.parse_args()
    file = args.file
    prefix = args.prefix
    count = args.count
    output_dir = os.path.abspath(args.output_dir)
    lines = args.lines

    # change logging format - https://realpython.com/python-logging/
    # format is - time, process_id, user, log level, message
    logging.basicConfig(level=args.loglevel, format='%(asctime)s - %(process)d - %(name)s - %(levelname)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')

    # check if file is not empty or get first 20 lines and see if it is a fasta file
    check_if_fasta(file,lines)

    # once confirmed fasta file, split the fasta
    split_fasta(file, prefix, count, output_dir)

def check_if_fasta(file,lines):
    """
    Check if input file is not empty and is in fasta format
    """
    # check if file is empty or not
    if os.stat(file).st_size == 0:
        logging.error(f"Looks like the file '{file}' is empty")
        sys.exit(1)

    fasta_format = False
    with open(file, 'r') as input_file:
        lines_cache = islice(input_file, lines)

        for current_line in lines_cache:
            current_line = current_line.rstrip("\n")
            if re.match(r'^>',current_line):
                fasta_format = True
    if not fasta_format:
        logging.error(f"""Looks like the input file is not in fasta format,
        please check the input file '{file}' or change the --lines parameter default '{lines}' value""")
        sys.exit(1)


def split_fasta(file, prefix, count, output_dir):
    """
    Split fasta file into chunks
    """
    # remove any chunks that already exists with same prefix
    remove_file(output_dir, prefix)
    chunk_counter = 1
    fasta_header_count = 1
    total_count = 0
    out_filename = prefix + "_" + str(chunk_counter) + ".txt"
    logging.info(f"Splitting fasta file into {count} fasta per chunk..")
    with open(file, 'r') as filehandle:
        for line in filehandle:
            line = line.rstrip("\n")
            # skip any blank lines
            if re.match(r'^\s*$',line):
                pass
            elif re.match(r'^>',line):
                total_count += 1
                # if chunk count is met refresh the chunk_counter
                if (fasta_header_count > count):
                    logging.debug(f"fasta counter met '{fasta_header_count}:{count}'")
                    # print(f"increment the counter from {chunk_counter} to {chunk_counter + 1}")
                    logging.debug(f"Created chunk - {out_filename}")
                    chunk_counter += 1
                    fasta_header_count = 1
                    out_filename = prefix + "_" + str(chunk_counter) + ".txt"
                    with open(out_filename, 'w') as out_file:
                        out_file.write(f"{line}\n")
                elif (fasta_header_count <= count):
                    # out_filename = prefix + "_" + str(chunk_counter) + ".txt"
                    with open(out_filename, 'a') as out_file:
                        out_file.write(f"{line}\n")
                        # break
                # increment the fasta counter
                fasta_header_count += 1
            else:
                with open(out_filename, 'a') as out_file:
                    out_file.write(f"{line}\n")
    logging.info(f"Total input fasta count:{total_count}")
    logging.info(f"Total chunks created:{int(total_count/count) + (total_count%count > 0)}") # https://stackoverflow.com/a/23590097
    logging.debug(f"Check: ({total_count}/{count}) = {int(total_count/count)}")
    logging.debug(f"Check: ({total_count}%{count} > 0) = {total_count%count > 0}")
    logging.debug(f"Total chunks created:{int(total_count/count) + total_count%count}, ({total_count} / {count})")




def remove_file(output_dir, prefix):
    """
    Remove files that are already present in the output directory
    """
    path_prefix = os.path.join(output_dir, prefix)
    for filename in glob.glob(path_prefix + "*"):
        try:
            logging.warn(f"Removing pre-exising chunks file before proceeding {filename}")
            os.remove(filename)
        except FileNotFoundError as err:
            logging.error(f"Error while deleting file {filename}. {err}")


if __name__ == "__main__":
    main()
