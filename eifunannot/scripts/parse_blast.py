#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to parse blast output

blast output generated in the format '-max_target_seqs 1 -evalue 1e-5 -outfmt \"6 qseqid sseqid pident qstart qend sstart send qlen slen length nident mismatch positive gapopen gaps evalue bitscore\"' is recommended


"""

# authorship and License information
__author__ = "Gemy George Kaithakottil"
__license__ = "GNU General Public License v3.0"
__maintainer__ = "Gemy George Kaithakottil"
__email__ = "gemygk@gmail.com"

import argparse
from argparse import RawTextHelpFormatter
from collections import defaultdict
import os
import re
import sys

# get script name
script = os.path.basename(sys.argv[0])

def merge_overlapping_intervals(coords):
    coords.sort(key=lambda interval: interval[0])
    merged = [coords[0]]
    for current in coords:
        previous = merged[-1]
        if current[0] <= previous[1]:
            previous[1] = max(previous[1], current[1])
        else:
            merged.append(current)
    return merged


def compute_blast_coverage(blast_tblr_output):
    # compute blast coverage
    blast_cov_info = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    blast_info = defaultdict(dict)
    for line in open(blast_tblr_output, encoding="utf8"):
        if line and not re.match(r"^\s*$", line) and not line.startswith("#"):
            line = line.rstrip("\n")
            x = line.split("\t")
            if len(x) != 17:
                raise ValueError(
                    "blast_tblr_output should have 17 columns, but provided {0}".format(
                        len(x)
                    )
                )
            (
                qseqid,
                sseqid,
                pident,
                qstart,
                qend,
                sstart,
                send,
                qlen,
                slen,
                length,
                nident,
                mismatch,
                positive,
                gapopen,
                gaps,
                evalue,
                bitscore,
            ) = x
            if int(qstart) > int(qend):
                qstart, qend = qend, qstart
            if int(sstart) > int(send):
                sstart, send = send, sstart
            # print(qseqid, qlen, qstart, qend, sseqid, slen, sstart, send)
            blast_info[qseqid]["qlen"] = qlen
            blast_info[qseqid]["slen"] = slen
            blast_cov_info[qseqid][sseqid]["qcoords"].append([int(qstart), int(qend)])
            blast_cov_info[qseqid][sseqid]["scoords"].append([int(sstart), int(send)])
    # print(blast_info)
    # print(blast_cov_info)

    # merge overlapping intervals
    for qseqid in blast_cov_info:
        for sseqid in blast_cov_info[qseqid]:
            qcoords = blast_cov_info[qseqid][sseqid]["qcoords"]
            scoords = blast_cov_info[qseqid][sseqid]["scoords"]
            qcoords = merge_overlapping_intervals(qcoords)
            scoords = merge_overlapping_intervals(scoords)
            blast_cov_info[qseqid][sseqid]["qcoords"] = qcoords
            blast_cov_info[qseqid][sseqid]["scoords"] = scoords

    # compute query and subject coverage
    for qseqid in blast_cov_info:
        qlen = int(blast_info[qseqid]["qlen"])
        for sseqid in blast_cov_info[qseqid]:
            slen = int(blast_info[qseqid]["slen"])
            qcoords = blast_cov_info[qseqid][sseqid]["qcoords"]
            scoords = blast_cov_info[qseqid][sseqid]["scoords"]
            qcov = 0
            for qcoord in qcoords:
                qcov += qcoord[1] - qcoord[0] + 1
            qper = f"{round(qcov / qlen * 100, 2):.2f}"
            scov = 0
            for scoord in scoords:
                scov += scoord[1] - scoord[0] + 1
            sper = f"{round(scov / slen * 100, 2):.2f}"
            blast_info[qseqid]["qcov"] = qcov
            blast_info[qseqid]["qper"] = qper
            blast_info[qseqid]["sseqid"] = sseqid
            blast_info[qseqid]["scov"] = scov
            blast_info[qseqid]["sper"] = sper
    return blast_info


def main():
    parser = argparse.ArgumentParser(
        description="Script to parse blast output",
        formatter_class=RawTextHelpFormatter,
        epilog="Note:\n"
        + "blast output generated in the format '-max_target_seqs 1 -evalue 1e-5 -outfmt \"6 qseqid sseqid pident qstart qend sstart send qlen slen length nident mismatch positive gapopen gaps evalue bitscore\"' is recommended"
        + "\n\nContact:"
        + __author__
        + "("
        + __email__
        + ")",
    )
    parser.add_argument(
        "-b",
        "--blast_tblr_output",
        help="Provide blast tabular output, see note below for recommended format",
        required=True,
    )
    args = parser.parse_args()
    blast_info = compute_blast_coverage(args.blast_tblr_output)
    # print header
    print("#qseqid", "#qlen", "#qcov", "#qcov_percent", "#sseqid", "#slen", "#scov", "#scov_percent", sep="\t")
    for qseqid in blast_info:
        print(
            qseqid,
            blast_info[qseqid]["qlen"],
            blast_info[qseqid]["qcov"],
            blast_info[qseqid]["qper"],
            blast_info[qseqid]["sseqid"],
            blast_info[qseqid]["slen"],
            blast_info[qseqid]["scov"],
            blast_info[qseqid]["sper"],
            sep="\t",
        )


if __name__ == "__main__":
    main()
