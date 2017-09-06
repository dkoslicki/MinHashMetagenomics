# Containment Min Hash
Containment min hash is a method that combines min hash and bloom filters that allows the Jaccard index (similarity of two sets defined as the cardinality of their intersection over their union) for sets of very different size.

This repository contains all the code necessary to re-construct [the associated paper](http://www.biorxiv.org/content/early/2017/09/04/184150).

## Requirements ##
This repository uses both Python2 code (for the [GemSIM](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3305602/) code) as well as Python3. See [the requirements](https://github.com/dkoslicki/MinHashMetagenomics/blob/master/src/Requirements.txt) document for which packages are needed (as well as other requirements).

## Usage ##
For the impatient, take a look at [MakePaper.sh](https://github.com/dkoslicki/MinHashMetagenomics/blob/master/src/MakePaper.sh) to see how the data is downloaded, the computations are run, and the paper is compiled.

The [main min hash implementation](https://github.com/dkoslicki/MinHashMetagenomics/blob/master/src/MinHash.py) is based off of [SourMash](https://github.com/dib-lab/sourmash) and is combined with [pybloom](https://github.com/jaybaird/python-bloomfilter) (or [pybloom_live](https://pypi.python.org/pypi/pybloom_live/) for Python3) for the simulated data, and combined with [jellyfish](https://github.com/gmarcais/Jellyfish) (a C implementation of bloom filters for DNA) for the real metagenomic data.
