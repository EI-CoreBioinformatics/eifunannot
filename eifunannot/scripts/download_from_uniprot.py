#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to download data from Swiss-Prot and Trembl databases in different formats

API Help:

https://www.ebi.ac.uk/training/events/programmatic-access-uniprot-using-python/
https://drive.google.com/file/d/1qZXLl19apibjszCXroC1Jg2BGjcychqX/view
https://colab.research.google.com/drive/1SU3j4VmXHYrYxOlrdK_NCBhajXb0iGes#scrollTo=fBpHl-UU1C3t

https://www.uniprot.org/help/api_queries  - USEFUL REFERENCE

https://www.uniprot.org/help/query-fields
https://www.uniprot.org/help/api
https://www.uniprot.org/help/api_retrieve_entries
https://www.uniprot.org/help/metalink
"""

# authorship and License information
__author__ = "Gemy George Kaithakottil"
__license__ = "GNU General Public License v3.0"
__maintainer__ = "Gemy George Kaithakottil"
__email__ = "gemygk@gmail.com"

# import libraries
import sys

try:
    assert sys.version_info.major >= 3
except AssertionError:
    sys.exit("Error: Python3 required, please 'source python_miniconda-4.8.3_py3.8_gk'")
import argparse
from argparse import RawTextHelpFormatter
import os
import re

import requests
from requests.adapters import HTTPAdapter, Retry
from datetime import datetime

# get script name
script = os.path.basename(sys.argv[0])

choices = ["fasta", "xlsx", "xml", "dat", "txt", "gff", "list"]

re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

date_fmt = datetime.now().strftime("%d_%m_%y_%H%M")


class DownloadFromUniprot:
    # Helper function to download data
    @staticmethod
    def get_next_link(headers):
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

    @staticmethod
    def get_batch_details(batch_url):
        while batch_url:
            response = session.get(batch_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            # print(f"headers:{response.headers}")
            release = response.headers["x-uniprot-release"]
            release_date = response.headers["x-uniprot-release-date"]
            download_date = response.headers["date"]
            return total, release, release_date, download_date

    @staticmethod
    def get_batch(batch_url):
        while batch_url:
            response = session.get(batch_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            yield response, total
            batch_url = DownloadFromUniprot.get_next_link(response.headers)

    def __init__(self, args):
        self.args = args
        self.taxon_name = args.taxon_name
        self.taxon_id = args.taxon_id
        self.format = args.format
        self.size = args.size
        self.progress = args.progress

    #  valid searches
    def download_from_uniprot(self, review_list, output):
        if self.format == "dat":
            url = f"https://rest.uniprot.org/uniprotkb/search?query=reviewed:{review_list}+AND+taxonomy_id:{self.taxon_id}&format=txt&size={self.size}"
        else:
            url = f"https://rest.uniprot.org/uniprotkb/search?query=reviewed:{review_list}+AND+taxonomy_id:{self.taxon_id}&format={self.format}&size={self.size}"
        progress = 0
        with open(output, "w") as f:
            for batch, total in DownloadFromUniprot.get_batch(url):
                # lines = batch.text.splitlines()
                lines = batch.text
                print(lines, file=f)

                if self.progress and self.format == "fasta":
                    count = 0
                    for line in lines.splitlines():
                        if line.startswith(">"):
                            count += 1
                    progress += count
                    print(f"{progress} / {total}")

    def run(self):
        for review_list in ["true", "false"]:
            if self.format == "dat":
                # +AND+taxonomy_name:{self.taxon_name}
                url = f"https://rest.uniprot.org/uniprotkb/search?query=reviewed:{review_list}+AND+taxonomy_id:{self.taxon_id}&format=txt&size={self.size}"
            else:
                url = f"https://rest.uniprot.org/uniprotkb/search?query=reviewed:{review_list}+AND+taxonomy_id:{self.taxon_id}&format={self.format}&size={self.size}"
            (
                total,
                release,
                release_date,
                download_date,
            ) = DownloadFromUniprot.get_batch_details(url)
            database = "SwissProt" if review_list == "true" else "TrEMBL"
            output = "_".join(
                [
                    "Uniprot",
                    database,
                    self.taxon_name,
                    self.taxon_id,
                    total,
                    date_fmt + "." + self.format,
                ]
            )
            print(
                f"# {download_date}\n"
                f"# Downloading {self.taxon_name}({self.taxon_id}) {database} {total} entries UniProt Release {release} (UniProt Release Date:{release_date}) to file '{output}' ... ",
                end="",
            )
            self.download_from_uniprot(review_list, output)
            print(f"done!")


def main():
    parser = argparse.ArgumentParser(
        description="Script to download data from Swiss-Prot and Trembl databases in different formats",
        formatter_class=RawTextHelpFormatter,
        epilog="where,"
        + "\n[format]        - Choose one of the [format] below:"
        + "\n                  fasta: fasta"
        # + "\n                  tab: Tab-seperated" #!TO DO
        + "\n                  xlsx: Excel" + "\n                  xml: XML"
        # + "\n                  rdf: RDF/XML" #!TO DO
        + "\n                  dat|txt: Text"
        + "\n                  gff: GFF"
        + "\n                  list: List"
        + "\n"
        + "\n    Note:"
        + "\n"
        + "\n    50557 is the taxon for class Insecta"
        + "\n    33090 is the taxon for class Viridiplantae"
        + "\n"
        + "\n    Example commands:"
        + "\n    download_from_uniprot Plasmodium_falciparum 5833 dat"
        + "\n    download_from_uniprot Plasmodium_falciparum 5833 txt"
        + "\n    download_from_uniprot Alveolata 33630 fasta"
        + "\n    download_from_uniprot Escherichia_coli 83333 fasta"
        + "\n\nContact:"
        + __author__
        + "("
        + __email__
        + ")",
    )
    parser.add_argument(
        "taxon_name",
        help="Give any name you like. Join the taxon name with underscore(s) if giving multiple words. Do NOT use quotes.",
    )
    parser.add_argument(
        "taxon_id",
        help="UniProt Taxon ID is the most importand field. Get it from UniProt first and use it here NOT from NCBI taxonomy ID",
    )
    parser.add_argument(
        "format",
        choices=choices,
        default="fasta",
        help="Choose one of the output file format (default: %(default)s)",
    )
    parser.add_argument(
        "--size",
        default=500,
        type=int,
        help="Always use size 500 as this will provide fast performance (default: %(default)s)",
    )
    parser.add_argument(
        "--progress",
        action="store_true",
        help="Report progress for downloads. Only available for fasta downloads (default: %(default)s)",
    )
    args = parser.parse_args()

    DownloadFromUniprot(args).run()


if __name__ == "__main__":
    main()
