# Containment Min Hash
Containment min hash is a method that combines min hash and bloom filters that allows the Jaccard index (similarity of two sets defined as the cardinality of their intersection over their union) for sets of very different size.

This repository contains all the code necessary to re-construct [the associated paper](http://www.biorxiv.org/content/early/2017/09/04/184150).

![alt text](https://github.com/dkoslicki/MinHashMetagenomics/blob/master/Paper/Figs/ClassicalConceptual.png "Cenceptual approach of min hash") ![alt text](https://github.com/dkoslicki/MinHashMetagenomics/blob/master/Paper/Figs/ContainmentConceptual.png "Cenceptual approach of containment min hash")


## Requirements ##
This repository uses both Python2 code (for the [GemSIM](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3305602/) code) as well as Python3. See [the requirements](https://github.com/dkoslicki/MinHashMetagenomics/blob/master/src/Requirements.txt) document for which packages are needed (as well as other requirements).

## Usage ##
For the impatient, take a look at [MakePaper.sh](https://github.com/dkoslicki/MinHashMetagenomics/blob/master/src/MakePaper.sh) to see how the data is downloaded, the computations are run, and the paper is compiled. You will need to change the variables at the top of [MakePaper.sh](https://github.com/dkoslicki/MinHashMetagenomics/blob/master/src/MakePaper.sh) to point to where you installed GemSIM, python2, etc.

The [main min hash implementation](https://github.com/dkoslicki/MinHashMetagenomics/blob/master/src/MinHash.py) is based off of [SourMash](https://github.com/dib-lab/sourmash) and is combined with [pybloom](https://github.com/jaybaird/python-bloomfilter) (or [pybloom_live](https://pypi.python.org/pypi/pybloom_live/) for Python3) for the simulated data, and combined with [jellyfish](https://github.com/gmarcais/Jellyfish) (a C implementation of bloom filters for DNA) for the real metagenomic data.

### Basic workflow ###
Here is a very basic example of how to implement the containment min hash approach:
```python
import MinHash as MH  # import the min hash package
from pybloom import BloomFilter  # import bloom filters for the larger set
import numpy as np

# Define variables
prime = 9999999999971  # taking hashes mod this prime
ksize = 11  # k-mer length
h = 10  # number of hashes to use
p = 0.01  # false positive rate for the bloom filter

# Create data
small_string = ''.join(np.random.choice(['A', 'C', 'T', 'G'], 100))  # small string to form the small set A
large_string = ''.join(np.random.choice(['A', 'C', 'T', 'G'], 1000)) + small_string  # large string to form the larger set B
A = set([small_string[i:i+ksize] for i in range(len(small_string) - ksize + 1)])  # the smaller set, in practice, you don't need to form this
B = set([large_string[i:i+ksize] for i in range(len(large_string) - ksize + 1)])  # the larger set, in practice, you don't need to form this
size_A = len(A)  # number of k-mers in A, this can be calculated when forming the sketch
size_B = len(B)  # number of k-mers in B, could use HyperLogLog on large_string to get the cardinality estimate instead

# Populate min hash sketch with smaller set
A_MH = MH.CountEstimator(n=h, max_prime=prime, ksize=ksize, save_kmers='y')
A_MH.add_sequence(small_string)  # creat the min hash of the small string

# Create the bloom filter and populate with the larger set
filt = BloomFilter(capacity=1.5*len(B), error_rate=p)  # Initialize the bloom filter
num_kmers = 0  #  used to count the number of k-mers in B
for kmer in B:
    if kmer not in filt:
        num_kmers += 1
        filt.add(kmer);

# Use the k-mers in the sketch of A and test if they are in the bloom filter of B
int_est = 0  # intersection estimate
for kmer in A_MH._kmers:
    if kmer is not '':  # in case the set was so small the Min Hash was not fully populated
        if kmer in filt:
            int_est += 1

int_est -= np.round(p*h)  # adjust for the false positive rate
containment_est = int_est / float(h)  # estimate of the containment index
jaccard_est = size_A * containment_est / (size_A + size_B - size_A * containment_est)
true_jaccard = len(A.intersection(B)) / float(len(A.union(B)))

print("Containment index estimate: %f" % containment_est)
print("Jaccard index estimate (via the containment approach): %f" % jaccard_est)
print("Relative error: %f" % (np.abs(jaccard_est-true_jaccard) / true_jaccard))
```

