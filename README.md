# BYASE

A library that uses Bayesian inference to identify gene-level 
and isoform-level ASE (Allele-specific expression) in polyploid 
(diploid or higher) organisms from single-end or 
paired-end RNA-seq data.

## Installation

BYASE is released as a Python package which requires **Python 3.6** 
or higher to be installed on the computer.

To install BYASE, some Python packages that BYASE depends on 
should be installed first. In order to to successfully compile 
some of the packages, several system libraries should be 
pre-installed. For example, on **Ubuntu 18.04**, these libraries 
may need to be installed:

```shell
sudo apt install zlib1g-dev libbz2-dev liblzma-dev
```

Then, use `pip` to install dependent packages:
```shell
pip3 install --user numpy scipy pandas cython
pip3 install --user pymc3 pyarrow
pip3 install --user pysam htseq
```

Then, use `pip` to install BYASE:
```shell
pip3 install --user byase
```

*After the installation, if you cannot run `byase` from the terminal, 
this is caused by the executable binary file `byase` not being found 
in the system path, you may need to run:*
```shell
export PATH=~/.local/bin:$PATH
``` 

To use the **plotting module**, `matplotlib` should also be 
installed:
```shell
pip3 install --user matplotlib
```
If BYASE is run on a server or cluster system, the plotting module 
may not be used, then the installation of package `matplotlib` can be 
ignored, because some problem may be encountered during the 
installation of `matplotlib` on such systems.


## Documentation

The documentation of BYASE can be found 
[here](https://byase-doc.readthedocs.io/en/latest/).
