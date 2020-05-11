#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Snakemake for a generic AHRD run
"""

# import modules
import os
import sys
import logging

# Request min version of snakemake
# https://snakemake.readthedocs.io/en/stable/snakefiles/writing_snakefiles.html#depend-on-a-minimum-snakemake-version
from snakemake.utils import min_version
min_version("5.9.1")

# declare variables
cwd = os.getcwd()

# get fasta
fasta = os.path.abspath(config["fasta"])
if not os.path.exists(fasta):
    print(f"ERROR: The fasta file cannot be accessed - '{fasta}'")
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

total_chunks = count / per_chunk
total_chunks = int(total_chunks) # avoid round-up

# check if there is remainder, then add one more to total chunks
if (count % per_chunk != 0):
    total_chunks += 1

print(f"INFO: Total number of fasta sequences:{count} [{fasta}]")
print(f"INFO: Total number of chunks:{total_chunks} [{per_chunk} per chunk]")
chunk_numbers = list(range(1,total_chunks+1)) # need to add chunks+1 to get desired length - check here https://stackoverflow.com/a/4504677

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
        expand(os.path.join(CHUNKS_FOLDER,"chunk_{sample}.txt"), sample=chunk_numbers),
        # blast database output
        expand(os.path.join(DATABASE_DIR,"{protein}.protein.fa.done"), protein=protein_samples),
        # blastp output
        expand(os.path.join(OUTPUT,"output_{protein}","chunk_{sample}.txt-vs-{protein}.blastp.{ext}"), sample=chunk_numbers, protein=protein_samples, ext=["tblr","completed"]),
        # interproscn output
        expand(os.path.join(INTERPROSCAN_DIR,"chunk_{sample}","chunk_{sample}.txt.interproscan.{ext}"), sample=chunk_numbers, ext=["tsv","completed"]),
        # ahrd output
        expand(os.path.join(AHRD_DIR,"chunk_{sample}","ahrd_output.{ext}"), sample=chunk_numbers, ext=["csv","completed"]),
        # collate ahrd output
        expand(os.path.join(OUTPUT,"ahrd_output.{ext}"), ext=["csv","completed"])

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
        basename = fasta_base
        # source = config["load"]["python"]
        # script = split_fasta_script
    shell:
        "(set +u" \
        + " && cd {params.cwd} " \
        # + " && {params.source} " \
        + " && ln -sf {input.fasta} "
        + " && /usr/bin/time -v split_fasta -v -f {params.basename} -p {params.prefix} -c {params.chunks}"
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
        output = os.path.join(OUTPUT,"output_{protein}","chunk_{sample}.txt-vs-{protein}.blastp.tblr"),
        completed = os.path.join(OUTPUT,"output_{protein}","chunk_{sample}.txt-vs-{protein}.blastp.completed")
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
        + " && /usr/bin/time -v blastp -db {input.database} -outfmt 6 -num_threads {params.threads} {params.parameters} -query {input.chunk} -out {output.output} " \
        + " && touch {output.completed} " \
        + ") 2> {log}"

# run interproscan_5_22_61
# -----------
rule interproscan_5_22_61:
    input:
        chunk = os.path.join(CHUNKS_FOLDER,"chunk_{sample}.txt"),
    output:
        output = os.path.join(INTERPROSCAN_DIR,"chunk_{sample}","chunk_{sample}.txt.interproscan.tsv"),
        completed = os.path.join(INTERPROSCAN_DIR,"chunk_{sample}","chunk_{sample}.txt.interproscan.completed"),
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
        # steps to clean input fasta
        + " && ln -sf {input} {params.temp_name} " \
        + " && {params.source_prinseq} " \
        + " && prinseq -fasta {params.temp_name} -aa -rm_header -out_good {params.temp_name}.good -out_bad {params.temp_name}.bad " \
        + " && mv {params.temp_name}.good.fasta {params.input} " \
        + " && {params.source_interproscan} " \
        + " && /usr/bin/time -v interproscan.sh -i {params.input} -b {params.prefix} -f TSV {params.parameters}" \
        + " && touch {output.completed}" \
        + ") 2> {log}"

# run ahrd
# -----------
rule ahrd:
    input:
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
        database_reference = os.path.join(DATABASE_DIR,"reference.protein.fa"),
        db_status_reference = os.path.join(DATABASE_DIR,"reference.protein.fa.done"),
        database_swissprot = os.path.join(DATABASE_DIR,"swissprot.protein.fa"),
        db_status_swissprot = os.path.join(DATABASE_DIR,"swissprot.protein.fa.done"),
        database_trembl = os.path.join(DATABASE_DIR,"trembl.protein.fa"),
        db_status_trembl = os.path.join(DATABASE_DIR,"trembl.protein.fa.done"),
        # blast results
        blast_reference = os.path.join(OUTPUT,"output_reference","chunk_{sample}.txt-vs-reference.blastp.tblr"),
        blast_swissprot = os.path.join(OUTPUT,"output_swissprot","chunk_{sample}.txt-vs-swissprot.blastp.tblr"),
        blast_trembl = os.path.join(OUTPUT,"output_trembl","chunk_{sample}.txt-vs-trembl.blastp.tblr"),
        # interproscan results
        ipr_results = os.path.join(INTERPROSCAN_DIR,"chunk_{sample}","chunk_{sample}.txt.interproscan.tsv"),
        # ahrd configuration
        ahrd_config = ahrd_config
    output:
        output = os.path.join(AHRD_DIR,"chunk_{sample}","ahrd_output.csv"),
        completed = os.path.join(AHRD_DIR,"chunk_{sample}","ahrd_output.completed")
    log:
        os.path.join(AHRD_DIR,"chunk_{sample}","ahrd.log")
    threads: 1
    params:
        cwd = os.path.join(AHRD_DIR,"chunk_{sample}"),
        source = config["load"]["ahrd"]
    shell:
        "(set +u" \
        + " && cd {params.cwd} " \
        + " && cp -a {input.blacklist_descline} {input.filter_descline_sprot}" \
        + " {input.filter_descline_trembl} {input.filter_descline_tair} {input.blacklist_token}" \
        + " {input.interpro_dtd} ." \
        + " && ln -sf {input.gene_ontology_result} {input.interpro_database} ." \
        + " && ln -sf {input.ahrd_config} ahrd_input_go_prediction.yml" \
        + " && ln -sf {input.chunk} proteins.fasta" \
        + " && ln -sf {input.ipr_results} interpro_result.raw" \
        + " && ln -sf {input.blast_reference} reference_blastp_tabular.txt" \
        + " && ln -sf {input.blast_swissprot} swissprot_blastp_tabular.txt" \
        + " && ln -sf {input.blast_trembl} trembl_blastp_tabular.txt" \
        + " && touch prepare_ahrd.done" \
        + " && {params.source} " \
        + " && /usr/bin/time -v java -Xmx10g -jar /ei/software/testing/ahrd/3.3.3/x86_64/bin/ahrd.jar ahrd_input_go_prediction.yml" \
        + " && touch {output.completed}" \
        + ") 2> {log}"


# run collate_ahrd
# -----------
rule collate_ahrd:
    input:
        expand(os.path.join(AHRD_DIR,"chunk_{sample}","ahrd_output.csv"),sample=chunk_numbers)
    output:
        output = os.path.join(OUTPUT,"ahrd_output.csv"),
        completed = os.path.join(OUTPUT,"ahrd_output.completed")
    log:
        os.path.join(AHRD_DIR,"collate_ahrd.log")
    threads: 1
    params:
        cwd = OUTPUT,
    shell:
        "(set +u" \
        + " && cd {params.cwd} " \
        + " && cat " + os.path.join(AHRD_DIR,"chunk_*","ahrd_output.csv") + " | awk '!/^#|^Protein-Accession|^$/' | sort -k1,1V > ahrd_output.woH.csv" \
        + " && head -n 3 " + os.path.join(AHRD_DIR,"chunk_1","ahrd_output.csv") + " | cat - ahrd_output.woH.csv > {output.output} " \
        + " && touch {output.completed}" \
        + ") 2> {log}"
