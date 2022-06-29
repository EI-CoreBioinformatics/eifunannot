#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to generate AHRD reference fasta from NCBI proteins

# Tuesday, 17 November 2020, 10:52AM

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
    """Get the header

    Returns:
        Fasta header field

    Examples:
        >>> get_header(">XP_003393040.1 bis(5'-nucleosyl)-tetraphosphatase [asymmetrical] [Bombus terrestris]")
        'XP_003393040.1'
        >>> get_header(">NP_001267818.1 acyl-CoA delta-9 desaturase [Bombus terrestris]")
        'NP_001267818.1'
    """
    header_id = re.search(r">([^\s]+)", line)  # get fasta header alone
    if header_id:
        header = header_id.group(1)
    else:
        logging.error(f"Cannot extract fasta header from line below:\n" f"{line}\n")
        sys.exit(1)
    return header


def get_description(line):
    """Get the description

    Returns:
        Description field

    Examples:
        >>> get_description(">XP_003393040.1 bis(5'-nucleosyl)-tetraphosphatase [asymmetrical] [Bombus terrestris]")
        "bis(5'-nucleosyl)-tetraphosphatase [asymmetrical]"
        >>> get_description(">NP_001267818.1 acyl-CoA delta-9 desaturase [Bombus terrestris]")
        'acyl-CoA delta-9 desaturase'
    """
    # >NP_001267818.1 acyl-CoA delta-9 desaturase [Bombus terrestris]
    # >XP_003393040.1 bis(5'-nucleosyl)-tetraphosphatase [asymmetrical] [Bombus terrestris]
    description_field = re.search(r">\S+\s+(.*)\s+\[", line)  # get description name
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
def format_fasta_header(fasta, keep_XP):
    """
    Format I need to make
    >ENSP00000479374.1 | Symbols:  | T cell-interacting, activating receptor on myeloid cells 1  |
    >ENSP00000487458.1 | Symbols:  | Unknown protein |
    >ENSP00000484596.1 | Symbols:  | Homo sapiens uncharacterized LOC389831 (LOC389831), transcript variant 3, mRNA.  |
    """

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
                # print curated protein with NP prefix
                if header.startswith("NP") and not keep_XP:
                    description = get_description(line)
                    if description != "Unknown function":
                        print_rest = True
                        print(
                            f">{header} | Symbols:  | {description.replace('[','(').replace(']',')')}  |"
                        )
                    else:
                        print_rest = False
                # print non-curated proteins
                elif keep_XP:
                    description = get_description(line)
                    if description != "Unknown function":
                        print_rest = True
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
        description="Script to generate AHRD reference fasta from NCBI proteins",
        formatter_class=RawTextHelpFormatter,
        epilog="Example command:\n"
        + script
        + " --fasta [ncbi.protein.faa]"
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
        help="Provide protein fasta file [ncbi.protein.faa]",
    )
    parser.add_argument(
        "--keep_XP",
        action="store_true",
        help="By default all XP_* fasta prefix are ignored and only NP_* fasta prefix protein descriptions are used (default: %(default)s)",
    )
    args = parser.parse_args()
    fasta = args.fasta
    keep_XP = args.keep_XP
    format_fasta_header(fasta, keep_XP)


if __name__ == "__main__":
    import doctest

    # doctest.testmod()
    # print(doctest.testmod())
    test_result = doctest.testmod()[0]
    # print(test_result)
    if test_result == 0:
        main()
    else:
        sys.exit()
    # https://stackoverflow.com/a/25691978
    # sys.exit(doctest.testmod()[0])

"""
# EXAMPLE
# GCF_000214255.1_Bter_1.0_protein.faa
>NP_001267818.1 acyl-CoA delta-9 desaturase [Bombus terrestris]
MAPNITSSPTGVLFEGETLDESSRVIDAPKTKYKRQIVWRNVIIFAYLHIGAVYGLYLALTSAKWATLLFALCLHLFSAL
GITAGAHRLWAHRSYKAKWPLQLFLMIGNTIAFQDAAIDWARDHRLHHKYSETNADPHNAKRGFFFAHVGWLLCRKHPDI
RAKGKGIDLSDLKNNPILSFQKKYYAILMPLLCFIVPTLIPVYCWDESWGNAYFVPTVLRYVYTLNMTWLVNSAAHMFGN
KPYDKYINPVENKMVAITALGEGWHNYHHVFPWDYKTAELGNYKVNITTLFIDACSKLGLAYDMKIVPQDLVRKRVERTG
DGSHNVWGWGDKDQTQQDRDVTMVVNLKKDH
>NP_001267823.1 venom serine protease precursor [Bombus terrestris]
MTGSKMLFACLALIAFLHPLVHVASAQECTTPDNKAGKCLGIRGCKPLLEMLQTQGHAAADFLRQSVCKYENNNPIVCCP
NEESREDRGILVEYEPLRPPHCGFSNVSHTRVVGGNPAVLGAWPWIAALGFRYPRNPALEPLWKCGGSLISSRHVLTAAH
CAEINELYVVRIGDLNLVRNDDGAHPVQIEIESKIIHPDYISGVTKHDIAILKLVEEVPFSEYVYPICLPVEDNLRNNNF
ERYYPFVAGWGSLAHHGPGSDDLMEVQVPVISNTECKNSYARFAAAHVTDTVLCAGYTQGGKDACQGDSGGPLMLPKKFT
FYQIGVVSYGHKCAAAGYPGVYTRVTSYLDDFILPAMQ


** NOTES ABOUT NP_ AND XP_ TAGS **

What is the difference between XM_ and NM_ accessions?
Accession numbers that begin with the prefix XM_ (mRNA), XR_ (non-coding RNA), and XP_ (protein) are model RefSeqs produced either by NCBI’s genome annotation pipeline or copied from computationally annotated submissions to the INSDC. These RefSeq records are derived from the genome sequence and have varying levels of transcript or protein homology support. They represent the predicted transcripts and proteins annotated on the NCBI RefSeq contigs and may differ from INSDC mRNA submissions or from the subsequently curated RefSeq records (with NM_, NR_, or NP_ accession prefixes). These differences may reflect real sequence variation (polymorphism), or errors or gaps in the available genome sequence. The support for model RefSeq records should be further evaluated by comparing them to other sequence information available in Gene, Related Sequences, and BLAST reports.
The genome annotation pipelines are automated and their predicted products may or may not be subject to manual curation, but the data may be refreshed periodically.
"""
