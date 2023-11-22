# Introduction

# Dependencies
  + bowtie2
  + primer3
  + python >=3.6:
    + primer3-py
    + biopython
    + beautifulsoup4
    + pysam
    + zipfile36

# installation
  ## Checkout repository
  ```
  #Grab repo from github
  git clone https://github.com/gofflab/HCRProbeDesign.git
  ```
  ## bowtie2 Installation
  To determine whether designed probesets target unique sequences within a given species, we use [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to quickly map all designed probesets against a reference genome index (mouse mm10 by default). Please ensure that bowtie2 is installed on your system to use this feature.  We recommend using conda to create a specific virtual environment and provide an `environment.yml` file to install bowtie2 and all other required dependencies for the probe design software.

  ```
  cd HCRProbeDesign/
  conda env create -f environment.yaml
  conda activate HCRProbeDesign
  ```

  ## pip install
  ```
  cd ..
  pip install ./HCRProbeDesign
  ```

# Quickstart
  ## Install mm10 bowtie2 index
  Before using `designProbes`, you must either disable genome masking `--no-genomemask` or install a bowtie2 index into the indices directory `<package_dir>/indices/`.
  For mouse (the default target species), we provide a simple command line utility `fetchIndex` to download and install the mouse mm10 bowtie2 directory into the package reference index folder.

  ```
  # Grab mm10 index and install in <package_install_dir>/indices/mm10/
  $ fetchIndex
  ```

# Command Line

`HCRProbeDesign` package provides the commandline executable `designProbes` as the primary
