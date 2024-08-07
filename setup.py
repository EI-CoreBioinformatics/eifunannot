# coding: utf-8
from setuptools import setup, find_packages
from setuptools.extension import Extension
from distutils.extension import Extension
from codecs import open
from os import path
import glob
import sys


here = path.abspath(path.dirname("__file__"))

with open(path.join(here, "DESCRIPTION.md"), encoding="utf-8") as description:
    description = description.read()

with open(path.join(here, "README.md"), encoding="utf-8") as fh:
    long_description = fh.read()

name = "eifunannot"
version = "1.5.1"

if sys.version_info.major != 3:
    raise EnvironmentError(
        """eifunannot is a python module that requires python3, and is not compatible with python2. Also, it is now 2020 and support for 2.x has ceased."""
    )


setup(
    name=name,
    version=version,
    author="Gemy Kaithakottil",
    author_email="Gemy.Kaithakottil@earlham.ac.uk",
    description=description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/EI-CoreBioinformatics/eifunannot",
    license="GPLv3",
    classifiers=[
        "Development Status :: 3 - Beta",
        "Topic :: Scientific Engineering :: Bio/Informatics",
        "License :: OSI Approved :: GPLv3 License",
        "Operating System :: POSIX :: Linux",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
    ],
    python_requires=">=3.6",
    zip_safe=False,
    keywords="gene functional annotation",
    # packages=find_packages(exclude=["test"]),
    packages=find_packages(),
    # scripts=glob.glob("eifunannot/scripts/split_*.py"),
    # scripts=[
    # 	script for script in glob.glob("bin/slurm/*_sub")
    # ],
    install_requires=[
        "snakemake>5.4.0",
        # "drmaa",
        "pandas",
    ],
    entry_points={
        "console_scripts": [
            "eifunannot=eifunannot.__main__:main",
            "split_fasta=eifunannot.scripts.split_fasta:main",
            "generate_ahrd_reference_fasta_from_ncbi=eifunannot.scripts.generate_ahrd_reference_fasta_from_ncbi:main",
            "generate_ahrd_reference_fasta_from_ensembl=eifunannot.scripts.generate_ahrd_reference_fasta_from_ensembl:main",
            "generate_ahrd_reference_fasta_from_file=eifunannot.scripts.generate_ahrd_reference_fasta_from_file:main",
            "download_from_uniprot=eifunannot.scripts.download_from_uniprot:main",
            "create_functional_annotation=eifunannot.scripts.create_functional_annotation:main",
            "parse_blast=eifunannot.scripts.parse_blast:main",
            "add_description_to_annotation_GFF3=eifunannot.scripts.add_description_to_annotation_GFF3:main",
        ]
    },
    package_data={
        "eifunannot.config": glob.glob("eifunannot/config/*json")
        + glob.glob("eifunannot/config/*yaml"),
        "eifunannot": ["*.smk"],
    },
    include_package_data=True,
)
