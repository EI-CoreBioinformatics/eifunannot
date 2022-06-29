#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to generate AHRD reference fasta from Ensembl proteins

# Tuesday, 12 May 2020, 12:55AM

"""


# import libraries
import argparse
from argparse import RawTextHelpFormatter
import os
import re
import sys
import logging

from eifunannot import __version__, __author__, __email__

# change logging format - https://realpython.com/python-logging/
# format is - time, process_id, user, log level, message
logging.basicConfig(
    format="%(asctime)s - %(process)d - %(name)s - %(levelname)s - %(message)s",
    datefmt="%d-%b-%y %H:%M:%S",
)

# check python version
try:
    assert sys.version_info >= (3, 7)
except AssertionError:
    logging.error(
        f"Python >=3.7 is required.\nCurrent is Python {sys.version}.\nExiting..."
    )
    sys.exit(1)

# get script name
script = os.path.basename(sys.argv[0])


def get_header(line):
    header_id = re.search(r">([^\s]+)", line)  # get fasta header alone
    if header_id:
        header = header_id.group(1)
    else:
        logging.error(f"Cannot extract fasta header from line below:\n" f"{line}\n")
        sys.exit(1)
    return header


def get_gene(line):
    # gene_id = re.search(r'>.*parent=(\w+),\w+',line)  # get gene name)
    gene_id = re.search(r">.*gene:([^\s]+)", line)  # get gene name
    if gene_id:
        gene = gene_id.group(1)
    else:
        logging.error(f"Cannot extract gene name from line below:\n" f"{line}\n")
        sys.exit(1)
    return gene


def get_gene_biotype(line):
    gene_biotype_field = re.search(
        r">.*gene_biotype:([^\s]+)", line
    )  # get gene_biotype name
    if gene_biotype_field:
        gene_biotype = gene_biotype_field.group(1)
    else:
        logging.error(f"Cannot extract gene_biotype from line below:\n" f"{line}\n")
        sys.exit(1)
    return gene_biotype


def get_transcript_biotype(line):
    transcript_biotype_field = re.search(
        r">.*transcript_biotype:([^\s]+)", line
    )  # get transcript_biotype name
    if transcript_biotype_field:
        transcript_biotype = transcript_biotype_field.group(1)
    else:
        logging.error(
            f"Cannot extract transcript_biotype from line below:\n" f"{line}\n"
        )
        sys.exit(1)
    return transcript_biotype


def get_description(line):
    description_field = re.search(
        r">.*description:([^\n]+)", line
    )  # get description name
    if description_field:
        description = description_field.group(1)
    else:
        logging.warning(
            f"Cannot extract description from protein_coding line below, assigning 'Unknown function'\n"
            f"{line}"
        )
        # sys.exit(1)
        description = "Unknown function"
    return description


# format_fasta_header
def format_fasta_header(fasta):
    """
    Format I need to make
    >ENSP00000479374.1 | Symbols:  | T cell-interacting, activating receptor on myeloid cells 1  |
    >ENSP00000487458.1 | Symbols:  | Unknown protein |
    >ENSP00000484596.1 | Symbols:  | Homo sapiens uncharacterized LOC389831 (LOC389831), transcript variant 3, mRNA.  |
    """

    # write the fasta information to file
    fasta_base = os.path.basename(fasta)
    # process the fasta
    print_rest = False
    with open(fasta, "r") as filehandle:
        for line in filehandle:
            line = line.rstrip("\n")
            # remove blank and comments
            if re.match(r"^\s*$", line) or line.startswith("#"):
                pass
            elif line.startswith(">"):
                header = get_header(line)
                gene = get_gene(line)
                gene_biotype = get_gene_biotype(line)
                transcript_biotype = get_transcript_biotype(line)
                # only worry about protein_coding at gene and transcript level
                if (
                    gene_biotype == "protein_coding"
                    and transcript_biotype == "protein_coding"
                ):
                    description = get_description(line)
                    if description != "Unknown function":
                        print_rest = True
                        # print("\t".join([header, gene, gene_biotype, transcript_biotype, description]))
                        # create format
                        # >ENSOABP00000000006.1 | Symbols:  | zgc:77880 (Source:ZFIN;Acc:ZDB-GENE-040426-1901)  |
                        # >ENSOABP00000000014.1 | Symbols:  | FERM domain containing 7 (Source:HGNC Symbol;Acc:HGNC:8079)  |
                        print(
                            f">{header} | Symbols:  | {description.replace('[','(').replace(']',')')}  |"
                        )
                    else:
                        print_rest = False
                else:
                    print_rest = False
            else:
                if print_rest:
                    print(line)


def main():
    parser = argparse.ArgumentParser(
        description="Script to generate AHRD reference fasta from Ensembl proteins",
        formatter_class=RawTextHelpFormatter,
        epilog="Example command:\n"
        + script
        + " --fasta [ensembl.pep.all.fa]"
        + "\n\nContact:"
        + __author__
        + "("
        + __email__
        + ")",
    )
    parser.add_argument(
        "-f",
        "--fasta",
        required=True,
        nargs="?",
        help="Provide protein fasta file [ensembl.pep.all.fa]",
    )
    args = parser.parse_args()
    fasta = args.fasta
    format_fasta_header(fasta)


if __name__ == "__main__":
    main()
