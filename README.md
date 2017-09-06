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
len_small_string = 100
len_large_string = 1000

# Create data: small set called A, large set called B
small_string = ''.join(np.random.choice(['A', 'C', 'T', 'G'], len_small_string))  # small string to form the small set A
size_A = len(set([small_string[i:i+ksize] for i in range(len(small_string) - ksize + 1)]))  # size of smaller set, used to convert containment index to Jaccard index
large_string = ''.join(np.random.choice(['A', 'C', 'T', 'G'], len_large_string)) + small_string  # large string to form the larger set B

# Populate min hash sketch with smaller set
A_MH = MH.CountEstimator(n=h, max_prime=prime, ksize=ksize, save_kmers='y')
A_MH.add_sequence(small_string)  # create the min hash of the small string

# Create the bloom filter and populate with the larger set
B_filt = BloomFilter(capacity=1.15*len_large_string, error_rate=p)  # Initialize the bloom filter
size_B_est = 0  # used to count the number of k-mers in B, could do much more intelligently (like with HyperLogLog)
for i in range(len(large_string) - ksize + 1):
	kmer = large_string[i:i+ksize]
	if kmer not in B_filt:
		size_B_est += 1
		B_filt.add(kmer)

# Use the k-mers in the sketch of A and test if they are in the bloom filter of B
int_est = 0  # intersection estimate
for kmer in A_MH._kmers:
	if kmer is not '':  # in case the set "A" was so small the Min Hash was not fully populated
		if kmer in B_filt:
			int_est += 1

int_est -= np.round(p*h)  # adjust for the false positive rate
containment_est = int_est / float(h)  # estimate of the containment index
jaccard_est = size_A * containment_est / (size_A + size_B_est - size_A * containment_est)

# calulate true jaccard for comparison
A = set([small_string[i:i+ksize] for i in range(len(small_string) - ksize + 1)])  # the smaller set
B = set([large_string[i:i+ksize] for i in range(len(large_string) - ksize + 1)])  # the larger set
size_A = len(A)  # number of k-mers in A
size_B = len(B)  # number of k-mers in B
true_jaccard = len(A.intersection(B)) / float(len(A.union(B)))

print("Containment index estimate: %f" % containment_est)
print("Jaccard index estimate (via the containment approach): %f" % jaccard_est)
print("True Jaccard index: %f" % true_jaccard)
print("Relative error: %f" % (np.abs(jaccard_est-true_jaccard) / true_jaccard))

```

