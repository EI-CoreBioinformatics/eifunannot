#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Snakemake for a generic AHRD run - rodents


Download the interproscan release specific datatabases
$ wget -N ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz
$ gunzip goa_uniprot_all.gaf.gz


- get the version of interproscan I have used
$ wget -N ftp://ftp.ebi.ac.uk/pub/databases/interpro/61.0/interpro.xml.gz
$ gunzip interpro.xml.gz

01/04/2019, 15:02:16
"""

# authorship and License information
__author__ = "Gemy George Kaithakottil"
__copyright__ = "Copyright 2018"
__license__ = "GNU General Public License v3.0"
__maintainer__ = "Gemy George Kaithakottil"
__email__ = "Gemy.Kaithakottil@gmail.com"
__status__ = "Production"
__version__ = "0.1"

# import modules
import os
import sys
import logging

# Request min version of snakemake
# https://snakemake.readthedocs.io/en/stable/snakefiles/writing_snakefiles.html#depend-on-a-minimum-snakemake-version
from snakemake.utils import min_version
min_version("5.4.0")

# declare variables
cwd = os.getcwd()

# get fasta
fasta = os.path.abspath(config["fasta"])
if not os.path.exists(fasta):
    print(f"ERROR: The fasta file cannot be accessed - '{fasta}'")
    # logging.error(f"ERROR: The fasta file cannot be accessed - '{fasta}'")
    sys.exit()
fasta_base = os.path.basename(fasta)

#######################
# FOLDERS
#######################
# get output folder
OUTPUT = os.path.abspath(config["output"])
# chunks folder
CHUNKS_FOLDER = os.path.join(OUTPUT,"data","chunks")
# blast database folder
DATABASE_DIR = os.path.join(OUTPUT,"data","database")
# interproscan folder
INTERPROSCAN_DIR = os.path.join(OUTPUT,"output_interproscan")
# ahrd output
AHRD_DIR = os.path.join(OUTPUT,"output_ahrd")
#######################

# AHRD config file
ahrd_config = os.path.abspath(config["ahrd_config"])

# need to add os.sep if we try to give absolute path to an already existing location
# here is the link explaining this issue
# https://stackoverflow.com/a/28080468
# res = os.path.join(os.sep,"ei","software","testing","ahrd","3.3.3","src","AHRD-3.3.3","test","resources")

# get proteins
protein_samples = []
all_protein_samples = config["databases"]
for protein_name, protein_path in all_protein_samples.items():
    protein_samples.append(protein_name)
    r_path = os.path.abspath(protein_path)
    if not os.path.exists(r_path):
        print (f"ERROR: '{protein_name}' fasta file '{protein_path}' cannot be accessed")
        # logging.error(f"ERROR: '{protein_name}' fasta file '{protein_path}' cannot be accessed")
        sys.exit()
    if not os.path.exists(DATABASE_DIR):
        os.makedirs(DATABASE_DIR)
    new_protein_name_list = (protein_name,"protein.fa")
    new_protein_name = ".".join(new_protein_name_list)
    cmd = "cd " + DATABASE_DIR + " && ln -sf " + r_path + " " + new_protein_name
    if not os.path.exists(os.path.join(DATABASE_DIR,new_protein_name)):
        os.system(cmd)

# create chunks
per_chunk = config["chunk_size"]
if not per_chunk:
    per_chunk = 500
    print(f'WARN: chunk_size option is required, use default values [{per_chunk}] instead')
    # logging.warning(f'WARN: chunk_size option is required, use default values [{per_chunk}] instead')
count = 0
with open(fasta, 'r') as file:
    for line in file:
        if line.startswith(">"):
            count += 1

# print("Total number of fasta sequences:{1} [{0}]".format(fasta,count))
print(f"INFO: Total number of fasta sequences:{count} [{fasta}]")
# logging.info(f"INFO: Total number of fasta sequences:{count} [{fasta}]")
total_chunks = count / per_chunk
total_chunks = int(total_chunks) # avoid round-up

# check if there is remainder, then add one more to total chunks
if (count % per_chunk != 0):
    total_chunks += 1

# print("Total number of chunks:{0} [{1} per chunk]".format(total_chunks,per_chunk))
print(f"INFO: Total number of chunks:{total_chunks} [{per_chunk} per chunk]")
# logging.info(f"INFO: Total number of chunks:{total_chunks} [{per_chunk} per chunk]")
# logging.warning("Total number of chunks:{0} [{1} per chunk]".format(total_chunks,per_chunk))
chunk_numbers = list(range(1,total_chunks+1)) # need to add chunks+1 to get desired length - check here https://stackoverflow.com/a/4504677
# print("Total chunk numbers:{}".format(chunk_numbers))
# FOR TESTING
# chunk_numbers = list(range(1,5)) # need to add chunks+1 to get desired length - check here https://stackoverflow.com/a/4504677


# summarise the inputs
# print ("Inputs provided:")
# print ("fasta:{0}".format(fasta))
print(f"INFO: Inputs provided:")
print(f"INFO: fasta:{fasta}")
# logging.info(f"INFO: Inputs provided:")
# logging.info(f"INFO: fasta:{fasta}")
# create logs folder
# need to find a proper fix for this as mentioned in the issue below, but for now using a quick fix
# # https://bitbucket.org/snakemake/snakemake/issues/838/how-to-create-output-folders-for-slurm-log#comment-45348663
cluster_logs_dir = os.path.join(cwd,"logs","cluster")
if not os.path.exists(cluster_logs_dir):
    os.makedirs(cluster_logs_dir)

#######################
# RULES STARTS HERE
#######################
shell.prefix("set -eo pipefail; ")

rule all:
    input:
        # chunk output
        expand(os.path.join(CHUNKS_FOLDER,"chunk_{sample}.txt"),sample=chunk_numbers),
        # blast database output
        expand(os.path.join(DATABASE_DIR,"{protein}.protein.fa.done"),protein=protein_samples),
        # blastp output
        expand(os.path.join(OUTPUT,"output_{protein}","chunk_{sample}.txt-vs-{protein}.blastp.tblr"),sample=chunk_numbers,protein=protein_samples),
        # interproscn output
        expand(os.path.join(INTERPROSCAN_DIR,"chunk_{sample}","chunk_{sample}.txt.interproscan.tsv"),sample=chunk_numbers),
        # ahrd output
        expand(os.path.join(AHRD_DIR,"chunk_{sample}","ahrd_output.csv"),sample=chunk_numbers),
        # collate ahrd output
        expand(os.path.join(OUTPUT,"ahrd_output.csv"))

#######################
# WORKFLOW
#######################

# run chunking
# ------------------
rule split_fasta:
    input:
        fasta = fasta
    output:
        expand(os.path.join(CHUNKS_FOLDER,"chunk_{sample}.txt"),sample=chunk_numbers)
    log:
        os.path.join(CHUNKS_FOLDER,"split_fasta.log")
    params:
        cwd = CHUNKS_FOLDER,
        prefix = "chunk",
        chunks = per_chunk,
        basename = fasta_base,
        source = config["load"]["python"]
    shell:
        "(set +u" \
        + " && cd {params.cwd} " \
        + " && {params.source} " \
        + " && ln -sf {input.fasta} "
        # + " && /usr/bin/time -v leaff_v0.2.pl {params.basename} {params.prefix} {params.chunks}"
        + " && /usr/bin/time -v /hpc-home/kaithakg/snakemake_scripts/AHRD_pipeline/0.1/split_fasta.v0.1.py -v -f {params.basename} -p {params.prefix} -c {params.chunks}"
        + ") 2> {log}"

# run blast makeblastdb
# -----------
rule makeblastdb:
    input:
        os.path.join(DATABASE_DIR,"{protein}.protein.fa")
    output:
        os.path.join(DATABASE_DIR,"{protein}.protein.fa.done")
    log:
        os.path.join(DATABASE_DIR,"makeblastdb.{protein}.log")
    threads: 1
    params:
        cwd = DATABASE_DIR,
        source = config["load"]["blast"]
    shell:
        "(set +u" \
        + " && cd {params.cwd} " \
        + " && {params.source} " \
        + " && /usr/bin/time -v makeblastdb -in {input} -dbtype prot && touch {output}" \
        + ") 2> {log}"

# run blastp
# -----------
rule blastp:
    input:
        chunk = os.path.join(CHUNKS_FOLDER,"chunk_{sample}.txt"),
        database = os.path.join(DATABASE_DIR,"{protein}.protein.fa"),
        db_status = os.path.join(DATABASE_DIR,"{protein}.protein.fa.done")
    output:
        os.path.join(OUTPUT,"output_{protein}","chunk_{sample}.txt-vs-{protein}.blastp.tblr")
    log:
        os.path.join(DATABASE_DIR,"blastp.chunk_{sample}_{protein}.log")
    params:
        cwd = OUTPUT,
        threads = "4",
        parameters = config["load_parameters"]["blast"],
        source = config["load"]["blast"]
    shell:
        "(set +u" \
        + " && cd {params.cwd} " \
        + " && {params.source} " \
        + " && /usr/bin/time -v blastp -db {input.database} -outfmt 6 -num_threads {params.threads} {params.parameters} -query {input.chunk} -out {output}" \
        + ") 2> {log}"

# run interproscan_5_22_61
# -----------
rule interproscan_5_22_61:
    input:
        chunk = os.path.join(CHUNKS_FOLDER,"chunk_{sample}.txt"),
    output:
        os.path.join(INTERPROSCAN_DIR,"chunk_{sample}","chunk_{sample}.txt.interproscan.tsv")
    log:
        os.path.join(DATABASE_DIR,"interproscan.chunk_{sample}.log")
    params:
        cwd = os.path.join(INTERPROSCAN_DIR,"chunk_{sample}"),
        temp_name = "chunk_{sample}.raw.txt",
        input = "chunk_{sample}.txt",
        prefix = "chunk_{sample}.txt.interproscan",
        threads = "4",
        parameters = config["load_parameters"]["interproscan"],
        source_prinseq = config["load"]["prinseq"],
        source_interproscan = config["load"]["interproscan"]
    shell:
        "(set +u" \
        + " && cd {params.cwd} " \
        # steps to remove asterix from the fasta
        + " && ln -sf {input} {params.temp_name} " \
        + " && {params.source_prinseq} " \
        + " && prinseq -fasta {params.temp_name} -aa -rm_header -out_good {params.temp_name}.good -out_bad {params.temp_name}.bad " \
        + " && mv {params.temp_name}.good.fasta {params.input} " \
        + " && {params.source_interproscan} " \
        # + " && /usr/bin/time -v interproscan.sh -i {params.input} -b chunk_{wildcards.sample}.txt.interproscan -f TSV {params.parameters}" \
        + " && /usr/bin/time -v interproscan.sh -i {params.input} -b {params.prefix} -f TSV {params.parameters}" \
        + ") 2> {log}"

# run ahrd
# -----------
rule ahrd:
    input:
        # lambda wildcards: config["databases"][wildcards.protein_new],
        chunk = os.path.join(CHUNKS_FOLDER,"chunk_{sample}.txt"),
        # from ahrd package
        blacklist_descline = os.path.abspath(config["blacklist_descline"]),
        filter_descline_sprot = os.path.abspath(config["filter_descline_sprot"]),
        filter_descline_trembl = os.path.abspath(config["filter_descline_trembl"]),
        filter_descline_tair = os.path.abspath(config["filter_descline_tair"]),
        blacklist_token = os.path.abspath(config["blacklist_token"]),
        interpro_dtd = os.path.abspath(config["interpro_dtd"]),
        # external data
        gene_ontology_result = os.path.abspath(config["gene_ontology_result"]),
        interpro_database = os.path.abspath(config["interpro_database"]),
        # other intputs from the earlier runs
        # database
        # database = os.path.join(DATABASE_DIR,"{protein}.protein.fa"),
        # db_status = os.path.join(DATABASE_DIR,"{protein}.protein.fa.done"),
        database_reference = os.path.join(DATABASE_DIR,"reference.protein.fa"),
        db_status_reference = os.path.join(DATABASE_DIR,"reference.protein.fa.done"),
        database_swissprot = os.path.join(DATABASE_DIR,"swissprot.protein.fa"),
        db_status_swissprot = os.path.join(DATABASE_DIR,"swissprot.protein.fa.done"),
        database_trembl = os.path.join(DATABASE_DIR,"trembl.protein.fa"),
        db_status_trembl = os.path.join(DATABASE_DIR,"trembl.protein.fa.done"),
        # blast results - NOT THE BEST WAY TO DO IT, BUT WILL NEED TO FIND A WAY TO DO IT PROPERLY
        # blast = os.path.join(OUTPUT,"output_{{protein}}","chunk_{sample}.txt-vs-{{protein}}.blastp.tblr"),
        blast_reference = os.path.join(OUTPUT,"output_reference","chunk_{sample}.txt-vs-reference.blastp.tblr"),
        blast_swissprot = os.path.join(OUTPUT,"output_swissprot","chunk_{sample}.txt-vs-swissprot.blastp.tblr"),
        blast_trembl = os.path.join(OUTPUT,"output_trembl","chunk_{sample}.txt-vs-trembl.blastp.tblr"),
        # interproscan results
        ipr_results = os.path.join(INTERPROSCAN_DIR,"chunk_{sample}","chunk_{sample}.txt.interproscan.tsv"),
        # ahrd configuration
        ahrd_config = ahrd_config
    output:
        os.path.join(AHRD_DIR,"chunk_{sample}","ahrd_output.csv")
    log:
        os.path.join(AHRD_DIR,"chunk_{sample}","ahrd.log")
    threads: 1
    params:
        cwd = os.path.join(AHRD_DIR,"chunk_{sample}"),
        source = config["load"]["ahrd"]
    shell:
        "(set +u" \
        + " && cd {params.cwd} " \
        + " && cp -a {input.chunk} {input.blacklist_descline} {input.filter_descline_sprot}" \
        + " {input.filter_descline_trembl} {input.filter_descline_tair} {input.blacklist_token}" \
        + " {input.interpro_dtd} ." \
        + " && ln -sf {input.gene_ontology_result} {input.interpro_database} ." \
        + " && ln -sf {input.ahrd_config} ahrd_input_go_prediction.yml" \
        + " && ln -sf {input.chunk} proteins.fasta" \
        + " && ln -sf {input.ipr_results} interpro_result.raw" \
        # + " && ln -sf {input.blast} {wildcards.protein}_blastp_tabular.txt" \ # TO DO
        + " && ln -sf {input.blast_reference} reference_blastp_tabular.txt" \
        + " && ln -sf {input.blast_swissprot} swissprot_blastp_tabular.txt" \
        + " && ln -sf {input.blast_trembl} trembl_blastp_tabular.txt" \
        + " && touch prepare_ahrd.done" \
        + " && {params.source} " \
        + " && /usr/bin/time -v java -Xmx10g -jar /ei/software/testing/ahrd/3.3.3/x86_64/bin/ahrd.jar ahrd_input_go_prediction.yml" \
        + ") 2> {log}"


# run collate_ahrd
# -----------
rule collate_ahrd:
    input:
        expand(os.path.join(AHRD_DIR,"chunk_{sample}","ahrd_output.csv"),sample=chunk_numbers)
    output:
        os.path.join(OUTPUT,"ahrd_output.csv")
    log:
        os.path.join(AHRD_DIR,"collate_ahrd.log")
    threads: 1
    params:
        cwd = OUTPUT,
    shell:
        "(set +u" \
        + " && cd {params.cwd} " \
        + " && cat " + os.path.join(AHRD_DIR,"chunk_*","ahrd_output.csv") + " | awk '!/^#|^Protein-Accession|^$/' | sort -k1,1V > ahrd_output.woH.csv" \
        + " && head -n 3 " + os.path.join(AHRD_DIR,"chunk_1","ahrd_output.csv") + " | cat - ahrd_output.woH.csv > {output} " \
        + ") 2> {log}"
