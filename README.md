
# PORI Cbioportal

This python package uses the IPR and GraphKB PORI adaptors to create PORI reports from dumps
of cbioportal data.

## Getting Started

### Install (For developers)

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

## Deployment (Publishing)

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

## Deployment (Scripts)

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
