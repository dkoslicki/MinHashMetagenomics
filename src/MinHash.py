"""
An implementation of a MinHash bottom sketch, applied to k-mers in DNA.
"""
from __future__ import print_function
import khmer
import screed
import h5py
import numpy as np
import os
import tempfile
import multiprocessing
from multiprocessing import Pool
import re
# To Do:
# Make the formation of the count vectors parallel

#CE = SM.CountEstimator(n=10, max_prime=1e10, ksize=10, input_filename='/Users/dkoslicki/Dropbox/Repositories/MinHash/data/Genomes/G000279475.fna')

notACTG = re.compile('[^ACTG]')

class CountEstimator(object):
    """
    A simple bottom n-sketch MinHash implementation.
    n is the number of sketches to keep
    Still don't know what max_prime is...
    """

    def __init__(self, n=None, max_prime=1e12, ksize=None, input_file_name=None, save_kmers='n', hash_list=None):
        if n is None:
            raise Exception
        if ksize is None:
            raise Exception

        if ksize % 2 == 0:
            raise Exception("Due to an issue with khmer, only odd ksizes are allowed")

        self.ksize = ksize
        self.hash_list = hash_list

        # get a prime to use for hashing
        p = get_prime_lt_x(max_prime)
        self.p = p

        # initialize sketch to size n
        #self._mins = [float("inf")]*n
        self._mins = [p]*n

        # initialize the corresponding counts
        self._counts = [0]*n

        # initialize the list of kmers used, if appropriate
        if save_kmers == 'y':
            self._kmers = ['']*n
        else:
            self._kmers = None

        # Initialize file name (if appropriate)
        self.input_file_name = input_file_name
        if self.input_file_name:
            self.parse_file()

    def parse_file(self):
        """
        opens a file and populates the CountEstimator with it
        """
        for record in screed.open(self.input_file_name):
            self.add_sequence(record.sequence)

    def add(self, kmer):
        """
        Add kmer into sketch, keeping sketch sorted, update counts accordingly
        """
        _mins = self._mins
        _counts = self._counts
        _kmers = self._kmers

        h = khmer.hash_murmur3(kmer)
        h = h % self.p
        if self.hash_list:  # If I only want to include hashes that occur in hash_list
            if h not in self.hash_list:  # If the kmer isn't in the hash_list, then break
                return

        if h >= _mins[-1]:
            return

        for i, v in enumerate(_mins):
            if h < v:
                _mins.insert(i, h)
                _mins.pop()
                _counts.insert(i, 1)
                _counts.pop()
                if self._kmers:
                    _kmers.insert(i, np.string_(kmer))
                    _kmers.pop()
                return
            elif h == v:
                _counts[i] += 1
                return
            # else: h > v, keep on going.

        assert 0, "should never reach this"

    def add_sequence(self, seq):
        """
         Sanitize and add a sequence to the sketch.
        """
        # seq = seq.upper().replace('N', 'G')
        seq = notACTG.sub('G', seq.upper())
        for kmer in kmers(seq, self.ksize):
            self.add(kmer)

    def jaccard_count(self, other):
        """
        Jaccard index weighted by counts
        """
        truelen = len(self._mins)
        while truelen and self._mins[truelen - 1] == self.p:  # If the value of the hash is the prime p, it doesn't contribute to the length of the hash
            truelen -= 1
        if truelen == 0:
            raise ValueError

        #return self.common(other) / float(truelen)
        (total1, total2) = self.common_count(other)
        return (total2 / float(sum(other._counts)), total1 / float(sum(self._counts)))
        # The entries here are returned as (A_{CE1,CE2}, A_{CE2,CE1})

    def jaccard(self, other):
        """
        Jaccard index
        """
        truelen = len(self._mins)
        while truelen and self._mins[truelen - 1] == self.p:  # If the value of the hash is the prime p, it doesn't contribute to the length of the hash
            truelen -= 1
        if truelen == 0:
            raise ValueError

        return self.common(other) / float(truelen)
    #similarity = jaccard_count

    def common_count(self, other):
        """
        Calculate number of common k-mers between two sketches, weighted by their counts
        """
        if self.ksize != other.ksize:
            raise Exception("different k-mer sizes - cannot compare")
        if self.p != other.p:
            raise Exception("different primes - cannot compare")

        common1 = 0
        common2 = 0
        for (count1, count2) in _yield_count_overlaps(self._mins, other._mins, self._counts, other._counts):
            common1 += count1  # The unpopulated hashes have count 0, so we don't have to worry about that here
            common2 += count2
        return (common1, common2)

    def common(self, other):
        """
        Calculate number of common k-mers between two sketches.
        """
        if self.ksize != other.ksize:
            raise Exception("different k-mer sizes - cannot compare")
        if self.p != other.p:
            raise Exception("different primes - cannot compare")

        common = 0
        for val in _yield_overlaps(self._mins, other._mins):
            if val != self.p:  # Make sure not to include the un-populated hashes p
                common += 1
        return common

    def _truncate(self, n):
        self._mins = self._mins[:n]

    def export(self, export_file_name):
        """
        This function will export the CountEstimator using hdf5
        """
        fid = h5py.File(export_file_name, 'w')
        grp = fid.create_group("CountEstimator")
        mins_data = grp.create_dataset("mins", data=self._mins)
        counts_data = grp.create_dataset("counts", data=self._counts)
        if self._kmers:
            kmer_data = grp.create_dataset("kmers", data=self._kmers)

        grp.attrs['class'] = np.string_("CountEstimator")
        grp.attrs['filename'] = np.string_(self.input_file_name)
        grp.attrs['ksize'] = self.ksize
        grp.attrs['prime'] = self.p
        fid.close()

    def jaccard_count_vector(self, other_list):
        """
        Function that returns the Y vector of MetaPalette. That is, the vector where Y[i] = Jaccard_count(self, other_CEs[i]
        :param other_list: a list of count estimator classes
        :return: a numpy vector with the same basis as other_list giving the jaccard_count of self with other[i]
        """
        Y = np.zeros(len(other_list))
        for i in range(len(other_list)):
            Y[i] = other_list[i].jaccard_count(self)[0]  # Gotta make sure it's not [1] (one's the CKM vector, the other is the "coverage")

        return Y

    def jaccard_vector(self, other_list):
        """
        Function that returns the Y vector of Jaccard values. That is, the vector where Y[i] = Jaccard(self, other_CEs[i]
        :param other_list: a list of count estimator classes
        :return: a numpy vector with the same basis as other_list giving the jaccard of self with other[i]
        """
        Y = np.zeros(len(other_list))
        for i in range(len(other_list)):
            Y[i] = other_list[i].jaccard(self)

        return Y


def import_single_hdf5(file_name):
    """
    This function will read an HDF5 file and populate the CountEstimator class accordingly
    :param file_name: input file name of HDF5 file created by CountEstimator.export(file_name)
    :return: CountEstimator
    """
    fid = h5py.File(file_name, 'r')  # This automatically handles non-existent files for me
    grp = fid["CountEstimator"]
    file_name = grp.attrs['filename']
    ksize = grp.attrs['ksize']
    prime = grp.attrs['prime']
    mins = list(grp["mins"][:])
    counts = list(grp["counts"][:])
    CE = CountEstimator(n=len(mins), max_prime=1e12, ksize=ksize)
    CE.p = prime
    CE._mins = mins
    CE._counts = counts
    CE.input_file_name = file_name
    if "kmers" in grp:
        CE._kmers = list(grp["kmers"][:])
    else:
        CE._kmers = None

    fid.close()
    return CE


def import_multiple_hdf5(input_files_list):
    """
    Import a bunch of HDF5 Count Estimators from a given list of HDF5 files
    :param file_names: List of HDF5 file names of Count Estimators
    :return: list of Count Estimators
    """
    CEs = list()
    #for file_name in input_files_list:
    #    CEs.append(import_single_hdf5(file_name))
    pool = Pool(processes=multiprocessing.cpu_count())
    CEs = pool.imap(import_single_hdf5, input_files_list, chunksize=130)
    pool.close()

    return CEs


def export_multiple_hdf5(CEs, out_folder):
    """
    Exports a list of Count Estimators to a bunch of HDF5 files in a certain folder
    :param CEs: a list of Count Estimators
    :return: None
    """
    for CE in CEs:
        if CE.input_file_name == None:
            raise Exception("This function only works when count estimator were formed from files (i.e. CE.input_filename != None")

    for CE in CEs:
        CE.export(os.path.join(out_folder, os.path.basename(CE.input_file_name) + ".CE.h5"))

    return


def export_multiple_to_single_hdf5(CEs, export_file_name):
    """
    This will take a list of count estimators and export them to a single, large HDF5 file
    :param CEs: list of Count Estimators
    :param file_name: output file name
    :return: None
    """
    fid = h5py.File(export_file_name, 'w')
    grp = fid.create_group("CountEstimators")
    for CE in CEs:
        subgrp = grp.create_group(os.path.basename(CE.input_file_name))  # the key of a subgroup is the basename (not the whole file)
        mins_data = subgrp.create_dataset("mins", data=CE._mins)
        counts_data = subgrp.create_dataset("counts", data=CE._counts)
        if CE._kmers:
            kmer_data = subgrp.create_dataset("kmers", data=CE._kmers)

        subgrp.attrs['class'] = np.string_("CountEstimator")
        subgrp.attrs['filename'] = np.string_(CE.input_file_name)  # But keep the full file name on hand
        subgrp.attrs['ksize'] = CE.ksize
        subgrp.attrs['prime'] = CE.p

    fid.close()


def import_multiple_from_single_hdf5(file_name, import_list=None):
    """
    This function will import multiple count estimators stored in a single HDF5 file.
    :param file_name: file name for the single HDF5 file
    :param import_list: List of names of files to import
    :return: a list of Count Estimators
    """
    CEs = list()
    fid = h5py.File(file_name, 'r')
    if "CountEstimators" not in fid:
        raise Exception("This function imports a single HDF5 file containing multiple sketches."
                        " It appears you've used it on a file containing a single sketch."
                        "Try using import_single_hdf5 instead")

    grp = fid["CountEstimators"]
    if import_list:
        iterator = [os.path.basename(item) for item in import_list]
    else:
        iterator = grp.keys()

    for key in iterator:
        if key not in grp:
            raise Exception("The key " + key + " is not in " + file_name)

        subgrp = grp[key]
        file_name = subgrp.attrs['filename']
        ksize = subgrp.attrs['ksize']
        prime = subgrp.attrs['prime']
        mins = list(subgrp["mins"][:])
        counts = list(subgrp["counts"][:])
        CE = CountEstimator(n=len(mins), max_prime=1e12, ksize=ksize)
        CE.p = prime
        CE._mins = mins
        CE._counts = counts
        CE.input_file_name = file_name
        if "kmers" in subgrp:
            CE._kmers = list(subgrp["kmers"][:])
        else:
            CE._kmers = None

        CEs.append(CE)

    fid.close()
    return(CEs)


class CE_map(object):
    """
    Helper function for mapping CountEstimator class over multiple input_file arguments
    """
    def __init__(self, n, max_prime, ksize, save_kmers):
        self.n = n
        self.max_prime = max_prime
        self.ksize = ksize
        self.save_kmers = save_kmers

    def __call__(self, input_file):
        return CountEstimator(n=self.n, max_prime=self.max_prime, ksize=self.ksize, input_file_name=input_file, save_kmers=self.save_kmers)


def compute_multiple(n=None, max_prime=1e12, ksize=None, input_files_list=None, save_kmers='n', num_threads=multiprocessing.cpu_count()):
    """
    Batch compute Count Estimators from a given list of file names.
    :param n: number of hashes to keep
    :param max_prime:
    :param ksize: kmer size to use
    :param input_files_list: list of input genomes (fasta/fastq)
    :param save_kmers: flag if you want to save kmers or not ('y' or 'n')
    :return: a list of Count Estimators
    """
    if n is None:
        raise Exception
    if ksize is None:
        raise Exception
    if input_files_list is None:
        raise Exception

    pool = Pool(processes=num_threads)
    CEs = pool.map(CE_map(n, max_prime, ksize, save_kmers), input_files_list)
    pool.close()
    return CEs


def form_common_kmer_matrix(CEs):
    """
    Forms the common kmer matrix for input list of count estimators
    :param CEs: a list of Count Estimators
    :return: a numpy matrix A where A_{i,j} \approx \sum_{w\in SW_k(g_i) \cap SW_k{g_j}} \frac{occ_w(g_j)}{|g_j| - k + 1}
    """
    A = np.zeros((len(CEs), len(CEs)), dtype=np.float64)
    for i in range(len(CEs)):
        for j in range(i, len(CEs)):
            (Aij, Aji) = CEs[i].jaccard_count(CEs[j])
            A[i, j] = Aij
            A[j, i] = Aji

    return A

def form_jaccard_kmer_matrix(CEs):
    """
    Forms a matrix A with A_{i,j} = Jaccard(CEs[i], CEs[j])
    :param CEs: list of Count Estimators
    :return: numpy matrix
    """
    A = np.zeros((len(CEs), len(CEs)), dtype=np.float64)
    for i in range(len(CEs)):
        for j in range(i, len(CEs)):
            val = CEs[i].jaccard_count(CEs[j])
            A[i, j] = val
            A[j, i] = val

    return A





##########################################



class CompositionSketch(object):
    def __init__(self, n=None, max_prime=1e12, ksize=None, prefixsize=None, input_file_name=None, save_kmers='n'):
        """
        :param n: the number of kmer hashes to keep in the sketch. Must be >= 1
        :param max_prime:
        :param ksize: the kmer size to use. Must be >= 1
        :param prefixsize: the size of the prefix that determines which hash to use. A total of 4**prefixsize hashes will be created. Must be  >= 1
        :return:
        """
        if n is None:
            raise Exception
        if ksize is None:
            raise Exception
        if prefixsize is None:
            raise Exception

        self.prefixsize = prefixsize
        self.ksize = ksize
        self.threshold = 0.01

        # get a prime to use for hashing
        p = get_prime_lt_x(max_prime)
        self.p = p

        # initialize 4**prefixsize MinHash sketches
        self.sketches = []
        for i in range(4**prefixsize):
            self.sketches.append(CountEstimator(n=n, max_prime=self.p, ksize=ksize, save_kmers=save_kmers))

        self.input_file_name = input_file_name
        if self.input_file_name:
            self.parse_file()

    def parse_file(self):
        """
        opens a file and populates the CountEstimator with it
        """
        for record in screed.open(self.input_file_name):
            self.add_sequence(record.sequence)

    def add(self, kmer):
        # This will choose the estimator/hash list based on the prefix of the kmer, then add the sketch of the full kmer to the estimator/hash list.
        #idx = khmer.forward_hash(kmer, self.prefixsize)  # This hashes just the first prefixsize number of nucleotides in kmer
        idx = 0
        for letter in kmer[0:self.prefixsize]:
            if letter == 'A':
                pass
            if letter == 'C':
                idx <<= 2
                idx += 1
            if letter == 'T':
                idx <<= 2
                idx += 2
            if letter == 'G':
                idx <<= 2
                idx += 3

        E = self.sketches[idx]
        E.add(kmer)
        #hash = khmer.hash_murmur3(kmer)
        #E.add(hash)

    def add_sequence(self, seq):
        "Sanitize and add a sequence to the sketch."
        # seq = seq.upper().replace('N', 'G')
        seq = notACTG.sub('G', seq.upper())
        for kmer in kmers(seq, self.ksize):
            self.add(kmer)

    def jaccard(self, other):
        """
        estimate of J(self,other) = |self \Cap other| / |self \Cup other|
        In actuality an average of the jaccards of the constituent Sketches
        """
        total = 0.
        count = 0
        for n, v in enumerate(self.sketches):
            try:
                total += v.jaccard(other.sketches[n])
                count += 1
            except ValueError:
                pass
        return total / float(count)

    def jaccard_count(self, other):
        """
        estimate of \sum_{w\in SW_k(g_i) \cap SW_k{g_j}} \frac{occ_w(g_j)}{|g_j| - k + 1}
        In actuality an average of the jaccards of the constituent Sketches
        """
        total = [0.,0.]
        count = 0
        for n, v in enumerate(self.sketches):
            try:
                res = v.jaccard_count(other.sketches[n])
                total[0] += res[0]
                total[1] += res[1]
                count += 1
            except ValueError:
                pass
        return (total[0] / float(count), total[1] / float(count))

    def similarity(self, other):
        """
        Returns the fraction of sketches that have jaccard > threshold
        """
        matches = 0
        count = 0
        for n, v in enumerate(self.sketches):
            try:
                f = v.jaccard(other.sketches[n])
                count += 1
                if f > self.threshold:
                    matches += 1
            except ValueError:
                pass
        return matches / float(count)

    def count_similarity(self, other):
        """
        Returns the fraction of sketches that have count_jaccard > threshold
        """
        matches = 0
        count = 0
        for n, v in enumerate(self.sketches):
            try:
                f = v.count_jaccard(other.sketches[n])
                count += 1
                if f > self.threshold:
                    matches += 1
            except ValueError:
                pass
        return matches / float(count)


def _yield_count_overlaps(mins1, mins2, counts1, counts2):
    """
    Return (\sum_{i \in indicies(mins1\cap min2)} counts1[i], \sum_{i \in indicies(mins1\cap min2)} counts2[i])
    """
    i = 0
    j = 0
    try:
        while 1:
            while mins1[i] < mins2[j]:
                i += 1
            while mins1[i] > mins2[j]:
                j += 1
            if mins1[i] == mins2[j]:
                yield (counts1[i], counts2[j])
                i += 1
                j += 1
    except IndexError:
        return


def _yield_overlaps(x1, x2):
    """yield common hash values while iterating over two sorted lists of hashes
    Returns an iterable object
    """
    i = 0
    j = 0
    try:
        while 1:
            while x1[i] < x2[j]:
                i += 1
            while x1[i] > x2[j]:
                j += 1
            if x1[i] == x2[j]:
                yield x1[i]
                i += 1
                j += 1
    except IndexError:
        return

def kmers(seq, ksize):
    """yield all k-mers of len ksize from seq.
    Returns an iterable object
    """
    for i in range(len(seq) - ksize + 1):
        yield seq[i:i+ksize]


# taken from khmer 2.0; original author Jason Pell.
def is_prime(number):
    """Check if a number is prime."""
    if number < 2:
        return False
    if number == 2:
        return True
    if number % 2 == 0:
        return False
    for _ in range(3, int(number ** 0.5) + 1, 2):
        if number % _ == 0:
            return False
    return True


def get_prime_lt_x(target):
    """Backward-find a prime smaller than (or equal to) target.

    Step backwards until a prime number (other than 2) has been
    found.

    Arguments: target -- the number to step backwards from
    """
    if target == 1:
        return 1

    i = int(target)
    if i % 2 == 0:
        i -= 1
    while i > 0:
        if is_prime(i):
            return i
        i -= 2

    if i <= 0:
        raise RuntimeError("unable to find a prime number < %d" % (target))


def test_jaccard_1():
    E1 = CountEstimator(n=0, ksize=21)
    E2 = CountEstimator(n=0, ksize=21)

    E1._mins = [1, 2, 3, 4, 5]
    E2._mins = [1, 2, 3, 4, 6]
    assert E1.jaccard(E2) == 4/5.0
    assert E2.jaccard(E1) == 4/5.0


def test_jaccard_2_difflen():
    E1 = CountEstimator(n=0, ksize=21)
    E2 = CountEstimator(n=0, ksize=21)

    E1._mins = [1, 2, 3, 4, 5]
    E2._mins = [1, 2, 3, 4]
    assert E1.jaccard(E2) == 4/5.0
    assert E2.jaccard(E1) == 4/4.0


def test_yield_overlaps():
    x1 = [1, 3, 5]
    x2 = [2, 4, 6]
    assert len(list(_yield_overlaps(x1, x2))) == 0


def test_yield_overlaps_2():
    x1 = [1, 3, 5]
    x2 = [1, 2, 4, 6]
    assert len(list(_yield_overlaps(x1, x2))) == 1
    assert len(list(_yield_overlaps(x2, x1))) == 1


def test_yield_overlaps_3():
    x1 = [1, 3, 6]
    x2 = [1, 2, 6]
    assert len(list(_yield_overlaps(x1, x2))) == 2
    assert len(list(_yield_overlaps(x2, x1))) == 2

def test_CompositionSketch():
    CS1 = CompositionSketch(n=5, max_prime=1e10, ksize=11, prefixsize=1, save_kmers='y')
    CS2 = CompositionSketch(n=5, max_prime=1e10, ksize=11, prefixsize=1, save_kmers='y')
    sequence1 = "CATGTGCATGTAGATCGATGCATGCATCGATGCATGATCGATCX"
    sequence2 = "CATGTGCATGTAGATCGATGCATGCATCGATGCATGATCGATCGAAAAAAAAAAAAAAAAAAA"
    CS1.add_sequence(sequence1)
    CS2.add_sequence(sequence2)
    assert CS1.similarity(CS2) == 1.0
    assert CS1.jaccard(CS2) == 0.6
    assert CS1.jaccard_count(CS2) == (0.5076923076923078, 0.5750000000000001)


def test_CountEstimator():
    CE1 = CountEstimator(n=5, max_prime=1e10, ksize=1)
    CE2 = CountEstimator(n=5, max_prime=1e10, ksize=1)
    sequence1 = "AAAAAAAA"  # 100% of the 1-mers of this sequence show up in the other
    sequence2 = "AAAACCCCCCCC"  # 4/12ths of the 1-mers in this sequence show up in the other
    CE1.add_sequence(sequence1)
    CE2.add_sequence(sequence2)
    assert CE1.jaccard_count(CE2) == (4/12., 1.0)
    assert CE2.jaccard_count(CE1) == (1.0, 4/12.)
    assert CE1.jaccard(CE2) == 1.0  # all of the unique kmers in seq1 show up in seq2
    assert CE2.jaccard(CE1) == 0.5  # half of the unique kmers in seq2 show up in seq1


def test_import_export():
    CE1 = CountEstimator(n=5, max_prime=1e10, ksize=1)
    CE2 = CountEstimator(n=5, max_prime=1e10, ksize=1)
    sequence1 = "AAAA"
    sequence2 = "AAAACCCC"
    CE1.add_sequence(sequence1)
    CE2.add_sequence(sequence2)
    temp_file = tempfile.mktemp()  # Make temporary file
    CE1.export(temp_file)  # Export the CountEstimator to temp file
    CE_Import = import_single_hdf5(temp_file)  # Read in the data
    os.remove(temp_file)  # Remove the temporary file
    assert CE_Import.jaccard_count(CE2) == CE1.jaccard_count(CE2)  # Make sure the results of the import are the same as the original CountEstimator


def test_hash_list():
    CE1 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y')
    seq1='acgtagtctagtctacgtagtcgttgtattataaaatcgtcgtagctagtgctat'
    CE1.add_sequence(seq1)
    hash_list = {424517919, 660397082L}
    CE2 = CountEstimator(n=5, max_prime=1e10, ksize=3, hash_list=hash_list, save_kmers='y')
    CE2.add_sequence(seq1)
    assert CE1.jaccard(CE2) == 0.4
    assert CE1.jaccard_count(CE2) == (1.0, 2/7.)


def test_vector_formation():
    CE1 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y')
    CE2 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y')
    CE3 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y')
    seq1 = 'tacgactgatgcatgatcgaactgatgcactcgtgatgc'
    seq2 = 'tacgactgatgcatgatcgaactgatgcactcgtgatgc'
    seq3 = 'ttgatactcaatccgcatgcatgcatgacgatgcatgatgtacgactgatgcatgatcgaactgatgcactcgtgatgczxerqwewdfhg'
    CE1.add_sequence(seq1)
    CE2.add_sequence(seq2)
    CE3.add_sequence(seq3)
    Y = CE1.jaccard_count_vector([CE1, CE2, CE3])
    assert (Y == np.array([1.,1.,0.5625])).all()
    Y2 = CE1.jaccard_vector([CE1, CE2, CE3])
    assert (Y2 == np.array([1.,1.,0.4])).all()


def test_suite():
    """
    Runs all the test functions
    :return:
    """
    test_jaccard_1()
    test_jaccard_2_difflen()
    test_yield_overlaps()
    test_yield_overlaps_2()
    test_yield_overlaps_3()
    test_CompositionSketch()
    test_CountEstimator()
    test_import_export()
    test_hash_list()
    test_vector_formation()

