#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Script to create functional annotation file

"""

# authorship and License information
__author__ = "Gemy George Kaithakottil"
__license__ = "GNU General Public License v3.0"
__maintainer__ = "Gemy George Kaithakottil"
__email__ = "gemygk@gmail.com"


# import libraries
import argparse
from argparse import RawTextHelpFormatter
import os
import re
import sys
from collections import defaultdict, namedtuple
from eifunannot.scripts.parse_blast import compute_blast_coverage

# get script name
script = os.path.basename(sys.argv[0])


# for blast
QUERY, QLEN, QCOV, QPER, TARGET, TLEN, TCOV, TPER = range(8)
Blast = namedtuple("Blast", "query qlen qcov qper target tlen tcov tper")

# for blast reference
QUERY, SYMBOL, DESC = range(3)
BlastRef = namedtuple("BlastRef", "query symbol desc")

# for ahrd
# Protein-Accession	Blast-Hit-Accession	AHRD-Quality-Code	Human-Readable-Description	Interpro-ID (Description)	Gene-Ontology-Term
ID, BLAST_HIT, AHRD_QC, HRD, IPRID, GO_TERM = range(6)
Ahrd = namedtuple("Ahrd", "id blast_hit ahrd_qc hrd iprid go_term")


# for metrics
TRANS, GENE, CONFIDENCE, BIOTYPE = range(4)
Metrics = namedtuple("Metrics", "trans gene confidence biotype")


# process blast
def process_blast(blast_tblr_output, blast_info):
    info = compute_blast_coverage(blast_tblr_output)
    for tid in info:
        blast_info[tid] = Blast(
            tid,
            info[tid]["qlen"],
            info[tid]["qcov"],
            info[tid]["qper"],
            info[tid]["sseqid"],
            info[tid]["slen"],
            info[tid]["scov"],
            info[tid]["sper"],
        )
    return blast_info

# process reference blast
def process_ref_blast(blast_reference_details, blast_ref_info):
    for line in open(blast_reference_details, encoding="utf8"):
        if line and not re.match(r"^\s*$", line) and not line.startswith("#"):
            line = line.rstrip("\n")
            # print(line)
            blast_ref_results = process_ref_blast_line(line)
            blast_ref_info[blast_ref_results.query] = blast_ref_results
    return blast_ref_info


def process_ref_blast_line(line):
    # QUERY, SYMBOL, DESC = range(3) # for reference
    fields = line.split("\t")
    results = BlastRef(
        fields[QUERY],
        fields[SYMBOL],
        fields[DESC],
    )
    return results


# process ahrd
def process_ahrd(ahrd_parsed_output, ahrd_info):
    for line in open(ahrd_parsed_output, encoding="utf8"):
        if (
            line
            and not re.match(r"^\s*$", line)
            and not line.startswith("#")
            and not line.startswith("Protein-Accession")
        ):  # ignore ahrd headers
            line = line.rstrip(
                "\n"
            )  # only strip with "\n" otherwise any empty tabs will be removed causing downstream issues
            ahrd_results = process_ahrd_line(line)
            ahrd_info[ahrd_results.id] = ahrd_results  # store to dict with the id
    return ahrd_info


def process_ahrd_line(line):
    # ID, BLAST_HIT, AHRD_QC, HRD, IPRID, GO_TERM = range(6) # for reference
    fields = line.split("\t")
    # if unknown the below variables are empty, so assign as None
    if not fields[BLAST_HIT]:
        fields[BLAST_HIT] = None
    if not fields[AHRD_QC]:
        fields[AHRD_QC] = None
    if not fields[IPRID]:
        fields[IPRID] = None
    if not fields[GO_TERM]:
        fields[GO_TERM] = None
    ahrd_results = Ahrd(
        fields[ID],
        fields[BLAST_HIT],
        fields[AHRD_QC],
        fields[HRD],
        fields[IPRID],
        fields[GO_TERM],
    )
    return ahrd_results


# process metrics
def process_metrics(metrics_parsed_output, metrics_output_cols, metrics_info):
    for line in open(metrics_parsed_output, encoding="utf8"):
        if (
            line
            and not re.match(r"^\s*$", line)
            and not line.startswith("#")
            and not line.startswith("#trans")
            and not line.startswith("TID")
        ):  # ignore metrics headers
            line = line.rstrip(
                "\n"
            )  # only strip with "\n" otherwise any empty tabs will be removed causing downstream issues
            metrics_results = process_metrics_line(metrics_output_cols, line)
            metrics_info[
                metrics_results.trans
            ] = metrics_results  # store to dict with the id
    return metrics_info


def process_metrics_line(metrics_output_cols, line):
    # TRANS, GENE, CONFIDENCE, BIOTYPE = range(4) # for reference
    fields = line.split("\t")
    metrics_results = Metrics(
        fields[metrics_output_cols[0]],
        fields[metrics_output_cols[1]],
        fields[metrics_output_cols[3]],
        fields[metrics_output_cols[2]],
    )
    return metrics_results


def main():
    parser = argparse.ArgumentParser(
        description="Script to create functional annotation file",
        formatter_class=RawTextHelpFormatter,
        epilog="Example command:\n"
        + script
        + " --ahrd_output ahrd_output.csv --blast_tblr_output blastp.tblr.overall_cov.txt"
        + "\n\nContact:"
        + __author__
        + "("
        + __email__
        + ")",
    )
    parser.add_argument(
        "--ahrd_output",
        required=True,
        help="Provide ahrd output file 'ahrd_output.csv'",
    )
    parser.add_argument(
        "--metrics_output",
        required=True,
        help="Provide mikado release metrics file *final_table.tsv",
    )
    parser.add_argument(
        "--metrics_output_cols",
        default="1,2,15,16",
        help="Provide mikado release metrics file columns to be used for TID,GID,biotype,confidence. Use quotes (default: %(default)s)",
    )
    parser.add_argument(
        "--blast_tblr_output",
        help="Provide parsed blast tabular file. blastp output generated in the format '-max_target_seqs 1 -evalue 1e-5 -outfmt \"6 qseqid sseqid pident qstart qend sstart send qlen slen length nident mismatch positive gapopen gaps evalue bitscore\"' is recommended",
    )
    parser.add_argument(
        "--blast_ref_name",
        default="Reference",
        help="Provide the reference name used for blast, use short and meaningful name (default: %(default)s)",
    )
    parser.add_argument(
        "--blast_reference_details",
        help="Provide a TSV file with additional functional information. Expect a three column TSV file - 'blast_ref_id symbol description' format (default: %(default)s)",
    )
    parser.add_argument(
        "--use_metrics_output_id",
        action="store_true",
        help="Enable this to use transcripts from metrics output file as reference, so that transposable_element_gene and ncrna_gene becomes part of the output. Default is to use transcripts from ahrd output as reference. (default: %(default)s)",
    )
    args = parser.parse_args()

    ahrd_output = args.ahrd_output
    metrics_output = args.metrics_output
    metrics_output_cols = args.metrics_output_cols.split(",")
    blast_tblr_output = args.blast_tblr_output
    blast_ref_name = args.blast_ref_name
    blast_reference_details = args.blast_reference_details
    use_metrics_output_id = args.use_metrics_output_id

    if len(metrics_output_cols) != 4:
        raise ValueError(
            "metrics_output_cols should have 4 columns, but provided {0}".format(
                len(metrics_output_cols)
            )
        )
    metrics_output_cols = [int(i) - 1 for i in metrics_output_cols]

    blast_info = {}  # create blast dictionary
    blast_ref_info = {}  # create blast reference dictionary
    if blast_tblr_output:
        # process blast
        blast_info = process_blast(blast_tblr_output, blast_info)
    if blast_reference_details:
        # process blast reference details
        blast_ref_info = process_ref_blast(blast_reference_details, blast_ref_info)

    # process ahrd
    ahrd_info = {}  # create ahrd dictionary
    ahrd_info = process_ahrd(ahrd_output, ahrd_info)

    # process ahrd
    metrics_info = {}  # create metrics dictionary
    metrics_info = process_metrics(metrics_output, metrics_output_cols, metrics_info)

    # pull together all the information
    header = ["#Transcipt", "#Gene", "#Confidence", "#Biotype"]

    if blast_tblr_output:
        header.extend(
            [
                f"#{blast_ref_name}-Blast-Hit",
                f"#Transcript-Blast-Hit-Coverage-By-{blast_ref_name}",
                f"#{blast_ref_name}-Blast-Hit-Coverage-By-Transcipt",
            ]
        )
    if blast_reference_details:
        header.extend(
            [
                f"#{blast_ref_name}-Blast-Hit-Symbol",
                f"#{blast_ref_name}-Blast-Hit-Description",
            ]
        )
    header.extend(
        [
            "#AHRD-Blast-Hit-Accession",
            "#AHRD-Quality-Code",
            "#Human-Readable-Description",
            "#Interpro-ID (Description)",
            "#Gene-Ontology-Term",
        ]
    )
    print(*header, sep="\t")

    for trans_id in (
        sorted(metrics_info) if use_metrics_output_id else sorted(ahrd_info)
    ):
        # ahrd = ahrd_info[trans_id]
        # print (ahrd)
        metrics = metrics_info[trans_id]
        ahrd = {}
        blast = Blast(None, None, None, None, None, None, None, None)
        blast_ref = BlastRef(None, None, None)
        if trans_id not in ahrd_info:
            # raise ValueError("Transcript not found in AHRD output file - '{0}'".format(trans_id))
            # sys.stderr(trans_id, "not found in Ahrd, defining all as empty")
            print(
                f"Warning: {trans_id} is not in AHRD, defining all as empty. Most likely these are not having protein sequences (like ncRNA's)",
                file=sys.stderr,
            )
            ahrd_results = Ahrd(trans_id, None, None, None, None, None)
            ahrd_info[ahrd_results.id] = ahrd_results
            ahrd = ahrd_info[trans_id]
        else:
            ahrd = ahrd_info[trans_id]

        if blast_tblr_output:
            if trans_id not in blast_info:
                # print (trans_id, "not found in blast, defining all as empty")
                blast_results = Blast(
                    trans_id, None, None, None, None, None, None, None
                )
                blast_info[blast_results.query] = blast_results
                blast = blast_info[trans_id]
            else:
                blast = blast_info[trans_id]

        if blast_reference_details:
            if blast.target not in blast_ref_info:
                # print(blast.target, "not found in blast, defining all as empty")
                blast_ref_results = BlastRef(blast.target, None, None)
                blast_ref_info[blast_ref_results.query] = blast_ref_results
                blast_ref = blast_ref_info[blast.target]
            else:
                blast_ref = blast_ref_info[blast.target]

        # now print
        output_line = [trans_id, metrics.gene, metrics.confidence, metrics.biotype]
        if blast_tblr_output:
            output_line.extend(
                [
                    blast.target,
                    blast.qper,
                    blast.tper,
                ]
            )
        if blast_reference_details:
            output_line.extend(
                [
                    blast_ref.symbol,
                    blast_ref.desc,
                ]
            )
        output_line.extend(
            [ahrd.blast_hit, ahrd.ahrd_qc, ahrd.hrd, ahrd.iprid, ahrd.go_term]
        )
        print(*output_line, sep="\t")


if __name__ == "__main__":
    main()
