#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

Script to generate AHRD reference fasta from a fasta and an annotation tsv file

"""

# import libraries
import argparse
from argparse import RawTextHelpFormatter
import os
import sys
import logging
from collections import defaultdict

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


class GenerateAHRDReferenceFasta(object):
    def __init__(self, args):
        self.args = args
        self.id_info = defaultdict(dict)

    def parse_file(self):
        with open(self.args.annotation_tsv, "r") as fh:
            for num, line in enumerate(fh, 1):
                line = line.rstrip("\n")
                # remove comments
                if line.startswith("#"):
                    continue
                # make sure that we atleast have more than two columns:
                # 1. id column
                # 2. description file
                x = line.split("\t")
                if len(x) < 3:
                    continue

                # get the main columns
                header = x[self.args.id - 1]
                if not header:
                    logging.warning(
                        f"Fasta header '{header}' is not defined in line '{num}' column number '{self.args.id}' below:"
                        f"\n{line}\nPlease correct the '{self.args.annotation_tsv}' file."
                    )
                    continue
                description = x[self.args.description - 1]
                if description == None or description == "":
                    description = "Unknown"
                symbol = x[self.args.symbol - 1] if self.args.symbol else None

                if header not in self.id_info:
                    self.id_info[header]["header"] = header
                    self.id_info[header]["description"] = description
                    self.id_info[header]["symbol"] = symbol
                else:
                    logging.error(
                        f"Duplicate header '{header}' encountered in line '{num}' below:"
                        f"\n{line}\nPlease correct the '{self.args.annotation_tsv}' file."
                    )
                    sys.exit(1)

    def format_fasta_header(self):
        """
        Format I need to make
        >YP_009370001.1 | Symbols: psbA | photosystem II protein D1 |
        >YP_009370002.1 | Symbols: matK | maturase K |
        >YP_009370004.1 | Symbols: psbK | photosystem II protein K |
        >YP_009370005.1 | Symbols: psbI | photosystem II protein I |
        """

        # write the fasta information to file
        fasta_base = os.path.basename(self.args.fasta)
        fasta_txt = fasta_base + ".ahrd_format.info.txt"
        if os.path.exists(fasta_txt):
            os.remove(fasta_txt)
        with open(fasta_txt, "w") as writer:
            writer.write(f"#Protein\t#Function\t#Symbol\n")

            # process the fasta
            with open(self.args.fasta, "r") as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if line.startswith(">"):
                        # get fasta header alone
                        try:
                            header = line.split(" ")[0].replace(">", "")
                        except:
                            raise ValueError(
                                f"Cannot extract fasta header from the fasta file '{self.args.fasta}' line below:\n{line}\n"
                            )
                        # check the gene in the fasta_info dictionary
                        if header in self.id_info:
                            if self.args.symbol:
                                print(
                                    f">{header} | Symbols: {self.id_info[header]['symbol']} | {self.id_info[header]['description']} |"
                                )
                            else:
                                print(
                                    f">{header} | Symbols:  | {self.id_info[header]['description']} ({self.id_info[header]['symbol']})  |"
                                )

                            writer.write(
                                f"{header}\t{self.id_info[header]['description']}\t{self.id_info[header]['symbol']}\n"
                            )
                        else:
                            logging.error(
                                f"Fasta header '{header}' information not present in the '{self.args.annotation_tsv}' file"
                            )
                    else:
                        print(line)
        logging.warning(f"#Fasta information is written to file '{fasta_txt}'")

    def run(self):
        self.parse_file()
        self.format_fasta_header()


def main():
    parser = argparse.ArgumentParser(
        description="Script to generate AHRD reference fasta from a fasta and an annotation tsv file",
        formatter_class=RawTextHelpFormatter,
        epilog="\n\nContact:" + __author__ + " (" + __email__ + ")",
    )
    parser.add_argument("fasta", help="Provide fasta file")
    parser.add_argument(
        "-a", "--annotation_tsv", required=True, help="Provide annotation tsv file"
    )
    parser.add_argument(
        "-i",
        "--id",
        required=True,
        type=int,
        help="Provide the column number in the annotation tsv file that match the input fasta header id",
    )
    parser.add_argument(
        "-d",
        "--description",
        required=True,
        type=int,
        help="Provide the column number in the annotation tsv file that match the functional description, for e.g., column in the annotation tsv file that has 'photosystem II protein D1'",
    )
    parser.add_argument(
        "-s",
        "--symbol",
        type=int,
        help="Provide the column number in the annotation tsv file that match the functional description symbol, for e.g., column in the annotation tsv file that has 'psbA'",
    )
    args = parser.parse_args()

    GenerateAHRDReferenceFasta(args).run()


if __name__ == "__main__":
    main()
