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
import logging
from collections import defaultdict
from pathlib import Path
from ftplib import FTP

# get script name
script = os.path.basename(sys.argv[0])

choices = ["fasta", "xlsx", "xml", "dat", "txt", "gff", "list"]

re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

date_fmt = datetime.now().strftime("%d_%m_%y_%H%M")

FORMAT = "# %(asctime)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)


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
        self.exclude_taxon_id = args.exclude_taxon_id
        self.filter_database = args.filter_database
        if self.filter_database:
            logging.info(
                "Change download '--format' to 'dat' as '--filter_database' option is enabled"
            )
            self.format = "dat"
        self.size = args.size
        self.progress = args.progress
        self.isoform_hash = defaultdict()

    #  valid searches
    def download_from_uniprot(self, url, output):
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

    def download_curated_isoforms(self):
        self.ftp = FTP("ftp.uniprot.org")
        # ftp.quit()
        logging.info("Connecting to UniProt ... ")
        self.ftp.login()
        logging.info("Done")

        logging.info("Download manually curated isoform sequences ... ")
        varsplic_out = "uniprot_sprot_varsplic.fasta.gz"
        self.ftp.cwd("/pub/databases/uniprot/current_release/knowledgebase/complete")
        with open(varsplic_out, "wb") as fp:
            self.ftp.retrbinary("RETR uniprot_sprot_varsplic.fasta.gz", fp.write)
        logging.info("Done")

        self.ftp.quit()

    def load_isoform_hash(self, isoform_file):
        acc = ""
        with open(isoform_file, "r") as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    # from:
                    # >sp|Q9S9Z8-2|14311_ARATH Isoform 2 of 14-3-3-like protein GF14 omicron OS=Arabidopsis thaliana OX=3702 GN=GRF11
                    # get:
                    # acc = Q9S9Z8
                    # fln script only matches to 9 splice variants with regex
                    # (^>\w+\|(\w+)\-\d\|.+)
                    # I am modifying it to all splice variants
                    # (^>\w+\|(\w+)\-\d+\|.+)

                    # for example, the above regex matches all below
                    # >sp|P05067-9|A4_HUMAN Isoform L-APP752 of Amyloid-beta precursor protein OS=Homo sapiens OX=9606 GN=APP
                    # >sp|P05067-10|A4_HUMAN Isoform APP639 of Amyloid-beta precursor protein OS=Homo sapiens OX=9606 GN=APP -- SKIPPED BY FLN RUBY SCRIPT
                    # >sp|P05067-11|A4_HUMAN Isoform 11 of Amyloid-beta precursor protein OS=Homo sapiens OX=9606 GN=APP -- SKIPPED BY FLN RUBY SCRIPT

                    """
                    All the splice variants gets added under the same 'acc'
                    """
                    acc_match_region = re.search(r"^>\w+\|(\w+)\-\d+\|.+", line)
                    if acc_match_region:
                        acc = acc_match_region.group(1)
                        if acc not in self.isoform_hash:
                            self.isoform_hash[acc] = f"{line}\n"
                        else:
                            self.isoform_hash[acc] += f"\n{line}\n"
                    else:
                        raise ValueError(
                            f"Error: Could extract relevant information from fasta header line '{line}' of file '{isoform_file}'"
                        )
                else:
                    self.isoform_hash[acc] += line

    def filter_incomplete_seqs(self, input_file, output):
        # UniProtKB fragments with FT NON_CONS and FT NON_TER features.
        #
        #     * FT NON_TER: The residue at an extremity of the sequence is not the terminal residue. If applied to position 1, this signifies that the first position is not the N-terminus of the complete molecule. If applied to the last position, it means that this position is not the C-terminus of the complete molecule. There is no description field for this key. Examples of NON_TER key feature lines:
        #       FT NON_TER 1 1
        #       FT NON_TER 29 29
        #     * FT NON_CONS: Non-consecutive residues. Indicates that two residues in a sequence are not consecutive and that there are a number of unreported or missing residues between them. Example of a NON_CONS key feature line:
        #       FT NON_CONS 1683 1684
        #
        # NON_CONS fragments are not indicated as non-consecutive in InterPro and being non-consecutive the match to methods may be incorrect if the method spans the 'break'.

        newseq = False
        print_seq = True
        id_id = ""
        acc_id = ""
        reviewed_status = ""
        description = ""
        organism_name = ""
        organism_taxid = ""
        gene_name = ""
        protein_existence = ""
        sequence_version = ""
        seq = ""
        organelle = ""

        # Open the file with writing permission
        output_file = open(output, "w")

        with open(input_file, "r") as fh:
            for line in fh:
                line = line.rstrip()
                if newseq is False:
                    id_match_region = re.search(r"^ID\s+(\w+)\s+(\w+);", line)
                    if id_match_region:
                        id_id = id_match_region.group(1)
                        reviewed_status = id_match_region.group(2)
                        newseq = True
                        description = ""
                        organism_name = ""
                        organism_taxid = ""
                        gene_name = ""
                        protein_existence = ""
                        sequence_version = ""
                        seq = ""
                        print_seq = True
                        organelle = ""

                else:
                    ac_match_region = re.search(r"^AC\s+(\w+);", line)
                    de_match_region = re.search(r"^DE\s+(.+);", line)
                    os_match_region = re.search(r"^OS\s+([^\(]+)", line)
                    ox_match_region = re.search(r"^OX\s+\w+=(\d+);", line)
                    gn_match_region = re.search(r"^GN\s+\w+=([^;]+)", line)
                    pe_match_region = re.search(r"^PE\s+(\d+):", line)
                    sv_match_region = re.search(r"^DT\s+.*sequence version\s+(\d+).", line)
                    og_match_region = re.search(r"^OG\s+(.+)", line)
                    seq_match_region = re.search(r"^\s+([\w\s]+)", line)
                    if ac_match_region:
                        acc_id = ac_match_region.group(1)
                    if de_match_region:
                        if description == "":
                            description = (
                                de_match_region.group(1)
                                .replace("RecName: Full=", "")
                                .replace("SubName: Full=", "")
                                .split("{")[0].strip()
                            )
                        if "Flags: Fragment" in line:
                            # print(f"#{acc_id} #{line}")
                            print_seq = False
                    elif os_match_region:
                        if organism_name == "":
                            organism_name = os_match_region.group(1).strip().rstrip(".")
                    elif ox_match_region:
                        organism_taxid = ox_match_region.group(1)
                    elif gn_match_region:
                        if gene_name == "":
                            gene_name = gn_match_region.group(1).split("{")[0].strip()
                    elif pe_match_region:
                        protein_existence = pe_match_region.group(1)
                    elif sv_match_region:
                        sequence_version = sv_match_region.group(1)
                    elif og_match_region:
                        organelle = og_match_region.group(1)
                    elif re.match(r"^FT\s+NON_TER\s+", line):
                        print_seq = False
                        # print(f"#{acc_id}   NON_TER")
                    elif re.match(r"^FT\s+NON_CONS\s+(\d+)\s+", line):
                        print_seq = False
                        # print(f"#{acc_id}   NON_CONS")
                    elif seq_match_region:
                        seq += seq_match_region.group(1)
                    elif re.match(r"^\/\/", line):
                        seq = re.sub(r"[\n\t\s]*", "", seq)
                        if seq[:1].lower() != "m":
                            print_seq = False
                        newseq = False

                        if print_seq:
                            header = ""
                            if reviewed_status == "Reviewed":
                                header += "sp|"
                            else:
                                header += "tr|"
                            header += f"{acc_id}|{id_id} {description} "
                            header += f"OS={organism_name} "
                            header += f"OX={organism_taxid} "
                            header += f"GN={gene_name} "
                            header += f"PE={protein_existence} "
                            header += f"SV={sequence_version}"

                            line_width = 60
                            sequence = "\n".join([seq[i:i+line_width] for i in range (0, len(seq), line_width)])
                            output_file.write(
                                f">{header}\n{sequence}\n"
                            )
                            # add the splice variants
                            if acc_id in self.isoform_hash:
                                output_file.write(f"{self.isoform_hash[acc_id]}\n")

    def run(self):
        # Download manually curated isoform sequences
        if self.filter_database:
            self.download_curated_isoforms()

            # extract the files
            logging.info("Extract downloaded files ... ")
            os.system("gunzip -f *.gz")
            logging.info("Done")

            # populate self.isoform_hash
            self.load_isoform_hash("uniprot_sprot_varsplic.fasta")
        for review_list in ["true", "false"]:
            # Uniprot Advanced: Share: Generate URI for API
            # https://rest.uniprot.org/uniprotkb/search?compressed=true&format=fasta&query=((taxonomy_id:3701)+NOT+(taxonomy_id:3702))+AND+(reviewed:true)&size=500
            url = f"https://rest.uniprot.org/uniprotkb/search?query=(reviewed:{review_list})+AND+((taxonomy_id:{self.taxon_id})"

            if self.exclude_taxon_id:
                url += f"+NOT+(taxonomy_id:{self.exclude_taxon_id})"

            if self.format == "dat":
                url += ")&format=txt"
            else:
                url += f")&format={self.format}"

            url += f"&size={self.size}"

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
            logging.info(f"Download date: {download_date}")
            message = f"Downloading taxon id ({self.taxon_id})"
            if self.exclude_taxon_id:
                message += f", excluding taxon id ({self.exclude_taxon_id}),"
            message += f" from {database} ({total} entries) UniProt Release {release} (UniProt Release Date:{release_date}) to file '{output}' ... "

            logging.info(message)
            self.download_from_uniprot(url, output)
            logging.info("Done")


            if self.filter_database:
                output_filtered = output.replace(".dat", ".filtered.fasta")
                logging.info(f"Filter incomplete sequences from '{output}' and add manually curated isoform sequences ... ")
                self.filter_incomplete_seqs(output, output_filtered)
                logging.info("Done")

                logging.info(f"Output file:'{output_filtered}'")


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
        + "\n"
        + "\n    # Download all 'flowering plants' (3398) but exclude 'bread wheat' (4565), plus filter database:"
        + "\n    download_from_uniprot Magnoliopsida_NOT_Triticum_aestivum 3398 fasta --exclude_taxon_id 4565 --filter_database"
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
        "--exclude_taxon_id",
        help="Provide UniProt Taxon ID you want to exclude. Get it from UniProt first and use it here NOT from NCBI taxonomy ID",
    )
    parser.add_argument(
        "--filter_database",
        action="store_true",
        help="Filter the database to remove incomplete sequences. If enabled the format of download will be in 'dat' format and filtered output will be in 'fasta' format (default: %(default)s)",
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
