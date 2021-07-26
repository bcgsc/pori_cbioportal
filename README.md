# PORI cBioportal

![build](https://github.com/bcgsc/pori_cbioportal/workflows/test/badge.svg) [![PyPi](https://img.shields.io/pypi/v/pori_cbioportal.svg)](https://pypi.org/project/pori_cbioportal)


This repository is part of the [Platform for Oncogenomic Reporting and Interpretation (PORI)](https://bcgsc.github.io/pori).


This python package uses the IPR and GraphKB PORI adaptors to create PORI reports from dumps
of cbioportal data.

## Getting Started

### Install via pip

```bash
pip install pori_cbioportal
```

### Download Study Data

Study data should be downloaded from cbioportal, for example

```bash
wget https://cbioportal-datahub.s3.amazonaws.com/laml_tcga_pan_can_atlas_2018.tar.gz
tar -xvzf laml_tcga_pan_can_atlas_2018.tar.gz
```

The folder should have the variant and metadata files, for example

```text
laml_tcga_pan_can_atlas_2018
|-- data_clinical_patient.txt
|-- data_clinical_sample.txt
|-- data_CNA.txt
|-- data_fusions.txt
|-- data_log2CNA.txt
|-- data_mutations_extended.txt
`-- data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt
```

### Generate Reports

This is then used to generate individual reports for all patients included in the study.
Note to do this you will need access to both a GraphKB server for matching and an IPR
server for upload.

```bash
pori_cbioportal laml_tcga_pan_can_atlas_2018 \
    --study_id "LAML TCGA" \
    --password $PASSWORD \
    --ipr_url https://YOUR_IPR_API_HOST/api \
    --graphkb_url https://YOUR_GRAPHKB_API_HOST/api
```

The loader will expect default names for the files but this can be overwritten with the other command line arguments. See the help menu for more options

```bash
pori_cbioportal --help
```

## Getting Started (For developers)

### Install

clone this repository

```bash
git clone ssh://git@svn.bcgsc.ca:7999/dat/pori_cbioportal.git
cd pori_cbioportal
```

create a virtual environment

```bash
python3 -m venv venv
source venv/bin/activate
```

install the package and its development dependencies

```bash
pip install -e .[dev]
```

Run the tests

```bash
pytest tests
```

### Deployment (Publishing)

Install the deployment dependencies

```bash
pip install .[deploy]
```

Build the distribution files

```bash
python setup.py install sdist bdist_wheel
```

Upload the distibutions to the package server (-r defined in your pypirc)

```bash
twine upload -r bcgsc dist/*
```

### Deployment (Scripts)

A buildout config is included by default which will build all console scripts defined
in the package.

create a virtual environment and install buildout

```bash
python3 -m venv venv
source venv/bin/activate
pip install -U pip setuptools zc.buildout
```

run buildout

```bash
buildout
```

This will create a directory `bin` with the executable scripts
