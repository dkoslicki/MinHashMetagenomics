# The first thing I am going to do is make sure I can implement the MinHash Sketch myself

#Virtual env in ~/KhmerForMinHash/bin/python

import khmer


# Installing hdf5 on mac in virtualenv
# download http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.16.tar.gz
# unpack, configure, make, sudo make install
# export HDF5_DIR=/Users/dkoslicki/KhmerForMinHash/hdf5-1.8.16/hdf5; /Users/dkoslicki/KhmerForMinHash/bin/pip install h5py




###########################################
# This was originally included in Titus' code, but I've modified it to CountEstimators
class Estimators(object):
    """
    A simple bottom n-sketch MinHash implementation.
    n is the number of sketches to keep
    Still don't know what max_prime is...
    """

    def __init__(self, n=None, max_prime=1e10, ksize=None):
        if n is None:
            raise Exception
        if ksize is None:
            raise Exception

        self.ksize = ksize

        # get a prime to use for hashing
        p = get_prime_lt_x(max_prime)
        self.p = p

        # initialize sketch to size n
        self._mins = [p]*n

    def add(self, kmer):
        "Add kmer into sketch, keeping sketch sorted."
        _mins = self._mins
        h = khmer.hash_murmur3(kmer)
        h = h % self.p  # This gurantees that the hash doesn't return the prime p

        if h >= _mins[-1]:
            return

        for i, v in enumerate(_mins):
            if h < v:
                _mins.insert(i, h)
                _mins.pop()
                return
            elif h == v:
                return
            # else: h > v, keep on going.

        assert 0, "should never reach this"

    def add_sequence(self, seq):
        "Sanitize and add a sequence to the sketch."
        seq = seq.upper().replace('N', 'G')
        for kmer in kmers(seq, self.ksize):
            self.add(kmer)

    def jaccard(self, other):
        truelen = len(self._mins)
        while truelen and self._mins[truelen - 1] == self.p:  # If the value of the hash is the prime p, it doesn't contribute to the length of the hash
            truelen -= 1
        if truelen == 0:
            raise ValueError

        return self.common(other) / float(truelen)
    similarity = jaccard

    def common(self, other):
        "Calculate number of common k-mers between two sketches."
        if self.ksize != other.ksize:
            raise Exception("different k-mer sizes - cannot compare")
        if self.p != other.p:
            raise Exception("different primes - cannot compare")

        common = 0
        for val in _yield_overlaps(self._mins, other._mins):
            common += 1
        return common

    def _truncate(self, n):
        self._mins = self._mins[:n]
