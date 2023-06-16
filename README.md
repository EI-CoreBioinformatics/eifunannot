# eifunannot
eifunannot - EI AHRD Functional Annotation Pipeline

eifunannot is a wrapper around [Automated Assignment of Human Readable Descriptions (AHRD)](https://github.com/groupschoof/AHRD) designed for execution in an HPC cluster (SLURM) environment.

Gemy Kaithakottil, David Swarbreck

## Getting Started

### Prerequisites

```Console
snakemake>=5.4.0
blast v2.6.0
prinseq v0.20.3
interproscan v5.22.61
ahrd v3.3.3
```
 
### Installation

First obtain the source code using

```console
$ cd /path/to/eifunannot/${version}
$ git clone git@github.com:EI-CoreBioinformatics/eifunannot.git src
$ cd src
# checkout the specific tag you want to install
$ git checkout tags/v1.4.0
```

To install, simply use from your current pip environment:
```console
$ version=1.4.0 && python setup.py bdist_wheel && pip install --prefix=/path/to/eifunannot/${version}/x86_64 -U dist/*whl
```

Create the wrapper

Make sure that both PATH and PYTHONPATH environments are updated for the new version
```console
cd /ei/software/cb/bin
$ cat eifunannot-1.4.0
#!/bin/bash
tool="eifunannot/1.4.0"
location="/ei/software/cb"
echo "${tool} is sourced from ${location} location"

source python_miniconda-4.7.12_py3.7_gk

echo "Usage:"
echo "  eifunannot --help"

export PYTHONPATH=/ei/software/cb/eifunannot/1.4.0/x86_64/lib/python3.7/site-packages
export PATH=/ei/software/cb/eifunannot/1.4.0/x86_64/bin:$PATH
```

### Run pipeline

Details to run the pipeline is provided in the [wiki page](../../wiki)
