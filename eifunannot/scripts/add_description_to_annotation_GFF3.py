#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to add description to release GFF3

"""

# authorship and License information
__author__ = "Gemy George Kaithakottil"
__license__ = "GNU General Public License v3.0"
__maintainer__ = "Gemy George Kaithakottil"
__email__ = "gemygk@gmail.com"


# import libraries
import argparse
import os
import re
import logging
import sys
from collections import defaultdict
from requests.utils import quote

# get script name
script = os.path.basename(sys.argv[0])

# get the GFF3 attributes
SEQID, SOURCE, TYPE, START, END, SCORE, STRAND, PHASE, ATTRIBUTE = range(9)

logging.basicConfig(
    format="%(asctime)s - %(process)d - %(name)s - %(levelname)s - %(message)s",
    datefmt="%d-%b-%y %H:%M:%S",
    level=logging.DEBUG,
)


class AddDescriptionToAnnotation:
    def __init__(self, args):
        self.args = args
        self.gff_file = args.gff_file
        self.annot_output = args.annot_output
        self.annot_column = args.annot_column
        self.source = args.source
        self.browser_specific = args.browser_specific

        self.gene_info = defaultdict()
        self.annot_info = defaultdict()

    def _encode(self, text):
        return quote(text, safe="")

    def process_annot_output(self):
        index = self.annot_column - 1
        with open(self.annot_output, "r") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if re.match(r"^\s*$", line) or line.startswith("#"):
                    continue
                x = line.split("\t")
                if len(x) < self.annot_column:
                    continue
                trans_id = x[0]
                if trans_id not in self.annot_info:
                    # print(trans_id, x[index])
                    self.annot_info[trans_id] = (
                        self._encode(x[index]) if self.browser_specific else x[index]
                    )

    def process_gff(self):
        with open(self.gff_file, "r") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if re.match(r"^\s*$", line):
                    continue
                if line.startswith("#"):
                    print(line)
                else:
                    row = line.split("\t")
                    if len(row) != 9:
                        raise ValueError(
                            f"Error: Not a standard GFF3 9 column line, see below:\n{line}"
                        )
                    # change source
                    if self.source:
                        row[SOURCE] = self.source
                    if row[TYPE].strip().lower() in [
                        "gene",
                        "ncrna_gene",
                        "pseudogene",
                    ]:
                        gattrib = {
                            k.lower(): v
                            for k, v in (
                                item.split("=")
                                for item in row[8].strip(" ;").split(";")
                            )
                        }
                        if any(
                            map(
                                lambda x: x is None,
                                (
                                    gattrib.get("ID".lower()),
                                    gattrib.get("Name".lower()),
                                    gattrib.get("biotype"),
                                    gattrib.get("confidence"),
                                ),
                            )
                        ):
                            raise ValueError(
                                "Error: Cannot parse all variables (ID, Name, biotype, confidence). It is required generating correct output. Please check entry:\n{}\n".format(
                                    "\t".join(row)
                                )
                            )
                        if gattrib["id"] not in self.gene_info:
                            self.gene_info[gattrib["id"]] = gattrib["biotype"]

                        if self.browser_specific:
                            row[TYPE] = "gene"
                        print(*row, sep="\t")

                    elif row[TYPE].strip().lower() in [
                        "mrna",
                        "ncrna",
                        "pseudogenic_transcript",
                    ]:
                        attrib = {
                            k.lower(): v
                            for k, v in (
                                item.split("=")
                                for item in row[8].strip(" ;").split(";")
                            )
                        }
                        if any(
                            map(
                                lambda x: x is None,
                                (
                                    attrib.get("ID".lower()),
                                    attrib.get("Parent".lower()),
                                    attrib.get("Name".lower()),
                                    attrib.get("Note".lower()),
                                    attrib.get("confidence"),
                                    attrib.get("representative"),
                                ),
                            )
                        ):
                            raise ValueError(
                                "Error: Cannot parse all variables (ID, Parent, Name, Note, confidence, representative). It is required generating correct output. Please check entry:\n{}\n".format(
                                    "\t".join(row)
                                )
                            )
                        if attrib["id"] not in self.annot_info:
                            logging.warning(
                                f"Looks like the transcript {attrib['id']} is not present in the input annotation file {self.annot_output}. Assigning as 'Unknown function'."
                            )
                            self.annot_info[attrib["id"]] = "Unknown function"
                        description = self.annot_info[attrib["id"]]
                        if attrib["parent"] not in self.gene_info:
                            raise ValueError(
                                f"Error: Looks like the input file is not sorted. The Parent {attrib['parent']} is not encountered before."
                            )
                        biotype = self.gene_info[attrib["parent"]]
                        new_name = "|".join(
                            [
                                f"{attrib['id']}",
                                f"{attrib['parent']}",
                                f"{biotype}",
                                f"conf:{attrib['confidence']}",
                                f"rep:{attrib['representative']}",
                            ]
                        )
                        new_attrib = ";".join(
                            [
                                f"ID={attrib['id']}",
                                f"Parent={attrib['parent']}",
                                f"Name={new_name}",
                                f"description={description}",
                            ]
                        )
                        if self.browser_specific:
                            row[TYPE] = "mRNA"
                            print(*row[:8], new_attrib, sep="\t")
                        else:
                            row[-1] = (
                                row[-1].rstrip(" ;") + f";description={description}"
                            )
                            print(*row, sep="\t")
                    else:
                        if self.browser_specific:
                            if row[TYPE] == "pseudogenic_exon":
                                row[TYPE] = "exon"
                        print(*row, sep="\t")

    def run(self):
        logging.info(f"Processing input file '{self.annot_output}'")
        self.process_annot_output()
        logging.info(f"Processing input file '{self.gff_file}'")
        self.process_gff()
        logging.info("Analysis complete")


class HelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter
):
    pass


def main():
    parser = argparse.ArgumentParser(
        prog=script,
        # formatter_class=lambda prog: argparse.HelpFormatter(prog, width=80),
        formatter_class=HelpFormatter,
        description="""
        Script to add description to release GFF3
        """,
        epilog=f"Contact: {__author__} ({__email__})",
    )
    parser.add_argument(
        "--gff_file",
        required=True,
        help="Provide minos.annotation.gff3 file",
    )
    parser.add_argument(
        "--annot_output",
        required=True,
        help="Provide annotation output tsv file",
    )
    parser.add_argument(
        "--annot_column",
        required=True,
        type=int,
        help="Provide annotation output column number to be used as description",
    )
    parser.add_argument(
        "--source",
        help="[optional] Provide a new source for the GFF3 file",
    )
    parser.add_argument(
        "--browser_specific",
        action="store_true",
        help="Enable this to change all type to standard format to load to apollo (default: %(default)s)",
    )
    args = parser.parse_args()
    AddDescriptionToAnnotation(args).run()


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        # Python flushes standard streams on exit; redirect remaining output
        # to devnull to avoid another BrokenPipeError at shutdown
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(1)  # Python exits with error code 1 on EPIPE
