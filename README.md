# BYASE

A library that uses Bayesian inference to identify gene-level 
and isoform-level ASE (Allele-specific expression) in polyploid 
(diploid or higher) organisms from single-end or 
paired-end RNA-seq data.

## Installation

BYASE is released as a Python package and can be installed by:
```shell
pip3 install --user byase
```
BYASE will automatically install some package dependencies, 
such as `numpy`, `scipy`, `pandas`, `pymc3`, `htseq` and 
`pyarrow`. In order to avoid possible installation problems 
on the server or cluster system, BYASE does not automatically 
install the package `matplotlib`. 
To use the **plotting module**, `matplotlib` should be 
installed manually on the computer:
```shell
pip3 install --user matplotlib
```

## Documentation

The documentation of BYASE can be found 
[here](https://byase-doc.readthedocs.io/en/latest/).
