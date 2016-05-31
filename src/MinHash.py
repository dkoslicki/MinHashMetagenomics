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
from itertools import *
import collections
from blist import *  # note, the import functions import the _mins etc. as lists, and the CE class imports them as blists.
# This shouldn't cause an issue, but will lead to slow performance if a CE is imported, then additional things are added.
# I.e. If you import a CE, don't add new elements, or you might have a bad day (or at least a long one).
import bisect
import scipy.optimize
import ctypes
import warnings
import subprocess
import filecmp
import shutil
import traceback
import random
warnings.simplefilter("ignore", RuntimeWarning)

# To Do:
# Implement hash_murmur3 to leave the khmer package. Need to implement reverse complement myself, etc.
# After that point, can use amino acids
# SNAP paired or single reads
# Get DIAMOND implemented
# Make the snap streaming automatically chunk the index_dirs if there are too many (can get max command len with xargs --show-limits)
# Make the count vector over a shared array (just like the kmer matricies)
# Use sam tools to partition the reads into aligned and unaligned (be careful with mate pairs)
# Implement paired reads for snap-aligner



notACTG = re.compile('[^ACTG]')


def unwrap_count_vector(arg):
    """
    Helper function for parallelizing the count_vector
    :param arg:
    :param kwarg:
    :return:
    """
    return arg[0].jaccard_count(arg[1])


def unwrap_jaccard_vector(arg):
    """
    Helper function for parallelizing the jaccard_vector
    :param arg:
    :param kwarg:
    :return:
    """
    return arg[0].jaccard(arg[1])


class CountEstimator(object):
    """
    A simple bottom n-sketch MinHash implementation.
    n is the number of sketches to keep
    Still don't know what max_prime is...
    """

    def __init__(self, n=None, max_prime=9999999999971., ksize=None, input_file_name=None, save_kmers='n', hash_list=None):
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
        self._mins = blist([p]*n)

        # initialize the corresponding counts
        self._counts = blist([0]*n)

        # initialize the list of kmers used, if appropriate
        if save_kmers == 'y':
            self._kmers = blist(['']*n)
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

        i = bisect.bisect_left(_mins, h)  # find index to insert h
        if _mins[i] == h:  # if h in mins, increment counts
            _counts[i] += 1
            return
        else:  # otherwise insert h, initialize counts to 1, and insert kmer if necessary
            _mins.insert(i, h)
            _mins.pop()
            _counts.insert(i, 1)
            _counts.pop()
            if _kmers:
                _kmers.insert(i, np.string_(kmer))
                _kmers.pop()
            return

        assert 0, "should never reach this"

    def add_sequence(self, seq):
        """
         Sanitize and add a sequence to the sketch.
        """
        # seq = seq.upper().replace('N', 'G')
        seq = notACTG.sub('G', seq.upper())  # more intelligent sanatization?
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

    def count_vector(self, other_list):
        """
        Function that returns the Y vector of MetaPalette. That is, the vector where Y[i] = Jaccard_count(self, other_CEs[i]
        :param other_list: a list of count estimator classes
        :return: a numpy vector with the same basis as other_list giving the jaccard_count of self with other[i]
        """
        Y = np.zeros(len(other_list))

        pool = Pool(processes=multiprocessing.cpu_count())
        Y_tuple = pool.map(unwrap_count_vector, zip([self] * len(other_list), other_list))
        pool.terminate()
        for i in range(len(other_list)):
            Y[i] = Y_tuple[i][1]  # Gotta make sure it's not [1] (one's the CKM vector, the other is the "coverage")

        return Y

    def jaccard_vector(self, other_list):
        """
        Function that returns the Y vector of Jaccard values. That is, the vector where Y[i] = Jaccard(self, other_CEs[i]
        :param other_list: a list of count estimator classes
        :return: a numpy vector with the same basis as other_list giving the jaccard of self with other[i]
        """
        pool = Pool(processes=multiprocessing.cpu_count())
        Y = np.array(pool.map(unwrap_jaccard_vector, zip([self]*len(other_list), other_list)))
        pool.terminate()

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
    mins = grp["mins"][...]  # For some reason, slicing is WAY slower than using ... in this case.
    counts = grp["counts"][...]
    CE = CountEstimator(n=len(mins), max_prime=3, ksize=ksize)
    CE.p = prime
    CE._mins = mins
    CE._counts = counts
    CE.input_file_name = file_name
    if "kmers" in grp:
        CE._kmers = grp["kmers"][...]
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
    pool = Pool(processes=multiprocessing.cpu_count())
    CEs = pool.map(import_single_hdf5, input_files_list, chunksize=144)
    pool.terminate()

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
        mins = subgrp["mins"][...]
        counts = subgrp["counts"][...]
        CE = CountEstimator(n=len(mins), max_prime=3, ksize=ksize)
        CE.p = prime
        CE._mins = mins
        CE._counts = counts
        CE.input_file_name = file_name
        if "kmers" in subgrp:
            CE._kmers = subgrp["kmers"][...]
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


def compute_multiple(n=None, max_prime=9999999999971., ksize=None, input_files_list=None, save_kmers='n', num_threads=multiprocessing.cpu_count()):
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


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]


def jaccard_count(ij):
    """
    Clone of jaccard_count from the count_estimator class, just so I can use shared memory arrays
    :param ij: a tuple of indicies to use in the global shared_mins and shared_counts
    :return: entries of the CKM matrix
    """
    mins1 = shared_mins[ij[0]]
    mins2 = shared_mins[ij[1]]
    counts1 = shared_counts[ij[0]]
    counts2 = shared_counts[ij[1]]
    truelen = len(mins1)
    while truelen and mins1[truelen - 1] == p:  # If the value of the hash is the prime p, it doesn't contribute to the length of the hash
        truelen -= 1
    if truelen == 0:
        raise ValueError
    common1 = 0
    common2 = 0
    for (count1, count2) in _yield_count_overlaps(mins1, mins2, counts1, counts2):
        common1 += count1  # The unpopulated hashes have count 0, so we don't have to worry about that here
        common2 += count2
    return (common2 / float(sum(counts2)), common1 / float(sum(counts1)))


def form_jaccard_count_matrix(all_CEs):
    """
    Forms the jaccard count kmer matrix when given a list of count estimators
    :param all_CEs: a list of count estimators
    :return: a numpy array of the jaccard count matrix
    """
    A = np.zeros((len(all_CEs), len(all_CEs)), dtype=np.float64)
    # Can precompute all the indicies
    indicies = []
    for i in xrange(len(all_CEs)):
        for j in xrange(len(all_CEs)):
            indicies.append((i, j))

    shared_mins_base = multiprocessing.Array(ctypes.c_double, len(all_CEs)*len(all_CEs[0]._mins))
    global shared_mins
    shared_mins = np.ctypeslib.as_array(shared_mins_base.get_obj())
    shared_mins = shared_mins.reshape(len(all_CEs), len(all_CEs[0]._mins))
    shared_counts_base = multiprocessing.Array(ctypes.c_double, len(all_CEs)*len(all_CEs[0]._counts))
    global shared_counts
    shared_counts = np.ctypeslib.as_array(shared_counts_base.get_obj())
    shared_counts = shared_counts.reshape(len(all_CEs), len(all_CEs[0]._counts))
    global p
    p = all_CEs[0].p
    for i in range(len(all_CEs)):
        shared_mins[i] = all_CEs[i]._mins
        shared_counts[i] = all_CEs[i]._counts

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    chunk_size = np.floor(len(indicies)/float(multiprocessing.cpu_count()))
    if chunk_size < 1:
        chunk_size = 1
    res = pool.imap(jaccard_count, indicies, chunksize=chunk_size)
    for (i, j), val in zip(indicies, res):
        A[i, j] = val[0]
        A[j, i] = val[1]

    pool.terminate()
    return A


def jaccard(ij):
    """
    Clone of jaccard_count from the count_estimator class, just so I can use shared memory arrays
    :param ij: a tuple of indicies to use in the global shared_mins and shared_counts
    :return: entries of the CKM matrix
    """
    mins1 = shared_mins[ij[0]]
    mins2 = shared_mins[ij[1]]
    truelen = len(mins1)
    while truelen and mins1[truelen - 1] == p:  # If the value of the hash is the prime p, it doesn't contribute to the length of the hash
        truelen -= 1
    if truelen == 0:
        raise ValueError

    common = 0
    for val in _yield_overlaps(mins1, mins2):
        if val != p:  # Make sure not to include the un-populated hashes p
            common += 1
    return common/float(truelen)


def form_jaccard_matrix(all_CEs):
    """
    Forms the jaccard count kmer matrix when given a list of count estimators
    :param all_CEs: a list of count estimators
    :return: a numpy array of the jaccard count matrix
    """
    A = np.zeros((len(all_CEs), len(all_CEs)), dtype=np.float64)
    # Can precompute all the indicies
    indicies = []
    for i in xrange(len(all_CEs)):
        for j in xrange(len(all_CEs)):
            indicies.append((i, j))

    shared_mins_base = multiprocessing.Array(ctypes.c_double, len(all_CEs)*len(all_CEs[0]._mins))
    global shared_mins
    shared_mins = np.ctypeslib.as_array(shared_mins_base.get_obj())
    shared_mins = shared_mins.reshape(len(all_CEs), len(all_CEs[0]._mins))
    global p
    p = all_CEs[0].p
    for i in range(len(all_CEs)):
        shared_mins[i] = all_CEs[i]._mins

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    chunk_size = np.floor(len(indicies)/float(multiprocessing.cpu_count()))
    if chunk_size < 1:
        chunk_size = 1
    res = pool.imap(jaccard, indicies, chunksize=chunk_size)
    for (i, j), val in zip(indicies, res):
        A[i, j] = val
        A[j, i] = val

    pool.terminate()
    return A


def lsqnonneg(A, Y, eps, machine_eps=1e-07):
    """
    Solve Ax=y for x using nonnegative least squares, with a threshold parameter eps. This is to be used if the A matrix is pre-computed
    :param A: The input common kmer or jaccard matrix
    :param Y: The input kmer or jaccard count vector
    :param eps: cutoff value. Only consider row/columns of A and entries of Y above this value
    :param machine_eps: All values below this are treated as zeros. Defaults to 1e-07
    :return: a vector x that approximately solves A x = Y subject to x >= 0 while ignoring rows/columns with Y[i] < eps
    """
    indicies = Y >= eps
    Y_eps = Y[indicies]
    A_eps = A[indicies, :][:, indicies]
    x_eps = scipy.optimize.nnls(A_eps, Y_eps)[0]
    x = np.zeros(len(Y))
    # repopulate the indicies
    x[indicies] = x_eps
    # set anything below machine_eps to zero
    for i in range(len(x)):
        if x[i] < machine_eps:
            x[i] = 0
    return x


def jaccard_count_lsqnonneg(CEs, Y, eps, machine_eps=1e-07):
    """
    Solve Ax=y for x using nonnegative least squares, with a threshold parameter eps. Only form A for subset of CEs with Y[i] > eps
    :param CEs: Input list of count estimators
    :param Y: The input kmer or jaccard count vector
    :param eps: cutoff value. Only consider row/columns of A and entries of Y above this value
    :param machine_eps: All values below this are treated as zeros. Defaults to 1e-07
    :return: a vector x that approximately solves A x = Y subject to x >= 0 while ignoring rows/columns with Y[i] < eps
    """
    indicies = Y >= eps
    if not indicies.any():
        raise Exception("No entries of Y are above eps: %f. Please decrease eps." % eps)
    Y_eps = Y[indicies]
    CEs_eps = []
    for i in range(len(indicies)):
        if indicies[i] == True:
            CEs_eps.append(CEs[i])

    A_eps = form_jaccard_count_matrix(CEs_eps)
    x_eps = scipy.optimize.nnls(A_eps, Y_eps)[0]
    # only take the entries of x_eps above the epsilon
    x_eps = x_eps/float(sum(x_eps))
    indicies_sorted = x_eps.argsort()[::-1]
    indicies_above_eps = []
    for index in indicies_sorted:
        indicies_above_eps.append(index)
        if x_eps[index] < eps:  # Add at leat one more index below the threshold (to make sure I don't have an empty list of indexes)
            break
    A_eps_above_eps = A_eps[indicies_above_eps, :][:, indicies_above_eps]
    x = np.zeros(len(Y))
    # repopulate the indicies
    x[indicies_above_eps] = x_eps
    # set anything below machine_eps to zero
    for i in range(len(x)):
        if x[i] < machine_eps:
            x[i] = 0
    return (x, A_eps_above_eps, indicies_above_eps)


def jaccard_lsqnonneg(CEs, Y, eps, machine_eps=1e-07):
    """
    Solve Ax=y for x using nonnegative least squares, with a threshold parameter eps. Only form A for subset of CEs with Y[i] > eps
    :param CEs: Input list of count estimators
    :param Y: The input kmer or jaccard count vector
    :param eps: cutoff value. Only consider row/columns of A and entries of Y above this value
    :param machine_eps: All values below this are treated as zeros. Defaults to 1e-07
    :return: a vector x that approximately solves A x = Y subject to x >= 0 while ignoring rows/columns with Y[i] < eps
    """
    indicies = Y >= eps
    Y_eps = Y[indicies]
    CEs_eps = []
    for i in range(len(indicies)):
        if indicies[i] == True:
            CEs_eps.append(CEs[i])

    A_eps = form_jaccard_matrix(CEs_eps)
    x_eps = scipy.optimize.nnls(A_eps, Y_eps)[0]
    # only take the entries of x_eps above the epsilon
    x_eps = x_eps/float(sum(x_eps))
    indicies_sorted = x_eps.argsort()[::-1]
    indicies_above_eps = []
    for index in indicies_sorted:
        indicies_above_eps.append(index)
        if x_eps[index] < eps:  # Add at leat one more index below the threshold (to make sure I don't have an empty list of indexes)
            break
    A_eps_above_eps = A_eps[indicies_above_eps, :][:, indicies_above_eps]
    x = np.zeros(len(Y))
    # repopulate the indicies
    x[indicies_above_eps] = x_eps
    # set anything below machine_eps to zero
    for i in range(len(x)):
        if x[i] < machine_eps:
            x[i] = 0
    return (x, A_eps_above_eps, indicies_above_eps)

##########################################



class CompositionSketch(object):
    def __init__(self, n=None, max_prime=9999999999971., ksize=None, prefixsize=None, input_file_name=None, save_kmers='n'):
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


def cluster_matrix(A_eps, A_indicies, taxonomy, cluster_eps=.01):
    """
    This function clusters the indicies of A_eps such that for a given cluster, there is another element in that cluster
    with similarity (based on A_eps) >= cluster_eps for another element in that same cluster. For two elements of
    distinct clusters, their similarity (based on A_eps) < cluster_eps.
    :param A_eps: The jaccard or jaccard_count matrix containing the similarities
    :param A_indicies: The basis of the matrix A_eps (in terms of all the CEs)
    :param cluster_eps: The similarity threshold to cluster on
    :return: (a list of sets of indicies defining the clusters, LCAs of the clusters)
    """
    #A_indicies_numerical = np.where(A_indicies == True)[0]
    A_indicies_numerical = A_indicies
    # initialize the clusters
    clusters = []
    for A_index in range(len(A_indicies_numerical)):
        # Find nearby elements
        nearby = set(np.where(A_eps[A_index, :] >= cluster_eps)[0]) | set(np.where(A_eps[:, A_index] >= cluster_eps)[0])
        in_flag = False
        in_counter = 0
        in_indicies = []
        for i in range(len(clusters)):
            if nearby & clusters[i]:
                clusters[i].update(nearby)  # add the nearby indicies to the cluster
                in_counter += 1  # keep track if the nearby elements belong to more than one of the previously formed clusters
                in_indicies.append(i)  # which clusters nearby shares elements with
                in_flag = True  # tells if it forms a new cluster
        if not in_flag:  # if new cluster, then append to clusters
            clusters.append(set(nearby))
        if in_counter > 1:  # If it belongs to more than one cluster, merge them together
            merged_cluster = set()
            for in_index in in_indicies[::-1]:
                merged_cluster.update(clusters[in_index])
                del clusters[in_index]  # delete the old clusters (now merged)
            clusters.append(merged_cluster)  # append the newly merged clusters
    clusters_full_indicies = []
    for cluster in clusters:
        cluster_full_indicies = set()
        for item in cluster:
            cluster_full_indicies.add(A_indicies_numerical[item])
        clusters_full_indicies.append(cluster_full_indicies)
    # Check to make sure the clustering didn't go wrong
    if sum([len(item) for item in clusters_full_indicies]) != len(A_indicies_numerical):  # Check the correct number of indicies
        raise Exception("For some reason, the total number of indicies in the clusters doesn't equal the number of indicies you started with")
    if set([item for subset in clusters_full_indicies for item in subset]) != set(A_indicies_numerical):  # Make sure no indicies were missed or added
        raise Exception("For some reason, the indicies in all the clusters doesn't match the indicies you started with")
    return clusters_full_indicies, cluster_LCAs(clusters_full_indicies, taxonomy)


def cluster_LCAs(clusters, taxonomy):
    """
    This function returns the lowest common ancestor in each one of the input clusters
    :param clusters: input clusters
    :param taxonomy: input taxonomy
    :return: a list with the ith element being the lowest common ancestor of the ith cluster
    """
    LCAs = []
    for cluster in clusters:
        found_LCA = False
        if len(cluster) == 1:
            LCAs.append(taxonomy[list(cluster)[0]].split()[2].split('|')[-1])
            found_LCA = True
            continue
        cluster_taxonomy = []
        for index in cluster:
            cluster_taxonomy.append(taxonomy[index])
        for rank in range(7, -1, -1):
            rank_names = []
            dummy_name = 0
            for tax_path in cluster_taxonomy:
                split_taxonomy = tax_path.split()[2].split('|')
                if len(split_taxonomy) < rank + 1:
                    rank_names.append(str(dummy_name))
                else:
                    rank_names.append(split_taxonomy[rank])
            if len(set(rank_names)) == 1 and "0" not in rank_names:
                LCAs.append(rank_names[0])
                found_LCA = True
                break
        if not found_LCA:
            LCAs.append('sk__-1_microorganism')  # In case they don't even have the kingdom in common
    return LCAs


def _write_single_cluster(tup):
    """
    Helper function. Writes a single fast file consisting of all the sequences in input_file_names
    :param tup: input tuple (out_dir, LCA, cluster_index, input_file_names)
    :return: the name of the created file
    """
    out_dir = tup[0]
    LCA = tup[1]
    cluster_index = tup[2]
    input_file_names = tup[3]
    out_file_name = os.path.join(out_dir, LCA + "_" + str(cluster_index) + "_" + ".fa")  # put the cluster index in the name in case there are shared LCAs
    out_file = open(out_file_name, 'w')
    i = 0
    for file_name in input_file_names:
        for record in screed.open(file_name):
            out_file.write(">" + LCA + "_" + str(i))
            out_file.write("\n")
            out_file.write(record.sequence)
            out_file.write("\n")
            i += 1
    out_file.close()
    return out_file_name

def make_cluster_fastas(out_dir, LCAs, clusters, CEs, threads=multiprocessing.cpu_count()):
    """
    This function will write a single fasta file for each of the clusters
    :param out_dir: the output directory in which to write the fasta files
    :param LCAs: the least common ancestors (from cluster_LCAs())
    :param clusters: the clusters (from cluster_matrix())
    :param CEs: The list of count estimators
    :param threads: number of threads to use
    :return: a list of files created (to be used in build_reference())
    """
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    pool = multiprocessing.Pool(processes=threads)
    file_names = pool.map(_write_single_cluster, zip(repeat(out_dir), LCAs, range(len(LCAs)), [[CEs[i].input_file_name for i in cluster] for cluster in clusters]), chunksize=1)
    pool.close()
    #pool.terminate()
    #pool.join()
    return file_names



##########################################################################
# SNAP aligner
def build_reference(reference_file, output_dir, large_index=True, seed_size=20, threads=multiprocessing.cpu_count(), binary="snap-aligner"):
    """
    Wrapper for the SNAP aligner index building
    :param reference_file: file from which an index will be made
    :param output_dir: directory in which to put the index files
    :param large_index: makes the index about 30% bigger, but faster for quick/inaccurate alignements
    :param seed_size: for initial match to begin alignment
    :param threads: number of threads to use
    :param binary: location of the snap-aligner binary
    :return: exit code of SNAP
    """
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    FNULL = open(os.devnull, 'w')
    if large_index:
        cmd = binary + " index " + reference_file + " " + output_dir + " -large -s " + str(seed_size) + " -t" + str(threads)
    else:
        cmd = binary + " index " + reference_file + " " + output_dir + " -s " + str(seed_size) + " -t" + str(threads)
    exit_code = subprocess.call(cmd, shell=True,  stdout=FNULL, stderr=subprocess.STDOUT)
    return exit_code


def _build_reference_helper(reference_file_name, output_dir, large_index, seed_size, threads, binary):
    prefix, ext = os.path.splitext(os.path.basename(reference_file_name))
    reference_dir = os.path.join(output_dir, prefix)
    try:
        exit_code = build_reference(reference_file_name, reference_dir, large_index=large_index, seed_size=seed_size, threads=threads, binary=binary)
    except:
        print('%s: %s' % (reference_file_name, traceback.format_exc()))
    return reference_dir


def _build_reference_star(args):
    return _build_reference_helper(*args)


def build_references(reference_files, output_dir, large_index=True, seed_size=20, threads=5, binary="snap-aligner"):
    """
    Will build indexes for SNAP aligner for multiple references
    :param reference_files: input reference fasta files
    :param output_dir: directory where all the indexes will be put
    :param large_index: makes the index about 30% bigger, but faster for quick/inaccurate alignements
    :param seed_size: for initial match to begin alignment
    :param threads: number of threads to use
    :param binary: location of the snap-aligner binary
    :return: list of location of all the indexes
    """
    if threads > 5:
        raise Exception("Python multiprocessing doesn't play well with many subprocesses. Please decrease the number of threads to something 5 or smaller.")
    pool = multiprocessing.Pool(processes=threads)
    index_dirs = pool.map(_build_reference_star, izip(reference_files, repeat(output_dir), repeat(large_index), repeat(seed_size), repeat(threads), repeat(binary)), chunksize=1)
    pool.close()
    pool.join()
    pool.terminate()
    return index_dirs


def build_one_reference_from_many(reference_files, output_dir, large_index=True, seed_size=20, threads=multiprocessing.cpu_count(), binary="snap-aligner"):
    # Make directory if necessary
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    # Concatenate all the reference files
    with open(os.path.join(output_dir, "all.fa"), 'w') as outfile:
        for fname in reference_files:
            with open(fname,'r') as infile:
                for line in infile:
                    outfile.write(line)
                infile.close()
        outfile.close()
    # Build the index
    exit_code = build_reference(os.path.join(output_dir, "all.fa"), output_dir, large_index=large_index, seed_size=seed_size, threads=threads, binary=binary)
    # Remove the concatenated file
    os.remove(os.path.join(output_dir, "all.fa"))
    return


def align_reads(index_dir, sample_file, out_file, filt='aligned', threads=multiprocessing.cpu_count(), edit_distance=14, min_read_len=50, binary="snap-aligner"):  # NOTE: snap-aligner will take SAM and BAM as INPUT!
    """
    Wrapper for the SNAP aligner.
    :param index_dir: Directory in which the index was placed (from build_reference
    :param sample_file: The file to align (can be fasta/fastq/sam/bam
    :param out_file: the output alignments. Can be *.sam or *.bam
    :param filt: If the output should be the aligned results ('aligned'), just the unaligned results ('unaligned'), or all results ('all')
    :param threads: number of threads to use
    :param edit_distance: the maximum edit distance to allow for an alignment
    :param min_read_len: Specify the minimum read length to align, reads shorter than this (after clipping) stay unaligned.  This should be
       a good bit bigger than the seed length or you might get some questionable alignments.  Default 50
    :param binary: location of the snap-aligner binary
    :return: exit code of SNAP
    """
    FNULL = open(os.devnull, 'w')
    if filt == 'aligned':
        #cmd = binary + " single " + index_dir + " " + sample_file + " -o " + out_file + " -f -F a -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len)
        cmd = binary + " single " + index_dir + " " + sample_file + " -o " + out_file + " -F a -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len)
    elif filt == 'unaligned':
        #cmd = binary + " single " + index_dir + " " + sample_file + " -o " + out_file + " -f -F u -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len)
        cmd = binary + " single " + index_dir + " " + sample_file + " -o " + out_file + " -F u -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len)
    elif filt == 'all':
        #cmd = binary + " single " + index_dir + " " + sample_file + " -o " + out_file + " -f -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len)
        cmd = binary + " single " + index_dir + " " + sample_file + " -o " + out_file + " -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len)
    else:
        raise Exception("aligned must be 'aligned', 'unaligned', or 'all'")
    exit_code = subprocess.call(cmd, shell=True,  stdout=FNULL, stderr=subprocess.STDOUT)
    FNULL.close()
    return exit_code


def stream_align_single(index_dirs, sample_file, out_file, format="bam", filt='all', threads=multiprocessing.cpu_count(), edit_distance=14, min_read_len=50, snap_binary="snap-aligner", samtools_binary="/local/cluster/samtools/bin/samtools"):
    """
    This will take a directory of indexes and stream filter the sample reads through it
    :param index_dirs: list of directories of snap indexes
    :param sample_file: sample file to align
    :param out_file: .sam or .bam file of reads that did/did not align to anything in the index_dirs
    :param filt: filtering option (unaligned, aligned, or all)
    :param threads: number of threads to use for each instance of snap-align
    :param edit_distance: max edit distance to allow for an alignment
    :param min_read_len: Specify the minimum read length to align, reads shorter than this (after clipping) stay unaligned.
     This should be a good bit bigger than the seed length or you might get some questionable alignments.  Default 50
    :param snap_binary: location of the snap-aligner snap_binary
    :return: exit code of SNAP
    """
    # cmd = snap_binary + " single " + index_dir + " " + sample_file + " -o -" + format + " - -f -F a -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len) + " | tee " + os.path.join(out_dir, os.path.basename(sample_file) + "_" + index_dir.split(os.sep)[-1] + "_" + filt + ".bam")
    if format != "bam" and format != "sam":
        raise Exception("Invalid format choice %s. Format must be one of 'bam' or 'sam'" % format)
    out_dir = os.path.dirname(out_file)
    FNULL = open(os.devnull, 'w')
    out_message_file_name = os.path.join(out_dir,"stream_align_stdout.txt")
    out_message_file = open(out_message_file_name, 'w')
    big_cmd = ''
    i = 0
    for index_dir in index_dirs:
        if i == 0:
            if filt == 'aligned':
                cmd = snap_binary + " single " + index_dir + " " + sample_file + " -o -" + format + " - -f -F a -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len)
            elif filt == 'unaligned':
                cmd = snap_binary + " single " + index_dir + " " + sample_file + " -o -" + format + " - -f -F u -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len)
            elif filt == 'all':
                if format == "bam":
                    cmd = snap_binary + " single " + index_dir + " " + sample_file + " -o -" + format + " - -f -xf 1.1 -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len) + " | " + samtools_binary + " view -@ 5 -h -b -f4 -o - -U " + os.path.join(out_dir, os.path.basename(sample_file) + "_" + index_dir.split(os.sep)[-1] + "_aligned.bam")
                else:
                    cmd = snap_binary + " single " + index_dir + " " + sample_file + " -o -" + format + " - -f -xf 1.1 -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len) + " | " + samtools_binary + " view -h -f4 -o - -U " + os.path.join(out_dir, os.path.basename(sample_file) + "_" + index_dir.split(os.sep)[-1] + "_aligned.sam")
            else:
                raise Exception("aligned must be 'aligned', 'unaligned', or 'all'")
            big_cmd = " " + cmd
        elif not i == len(index_dirs)-1:
            if filt == 'aligned':
                cmd = snap_binary + " single " + index_dir + " -" + format + " - -o -" + format + " - -f -F a -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len)
            elif filt == 'unaligned':
                cmd = snap_binary + " single " + index_dir + " -" + format + " - -o -" + format + " - -f -F u -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len)
            elif filt == 'all':
                if format == "bam":
                    cmd = snap_binary + " single " + index_dir + " -" + format + " - -o -" + format + " - -f -xf 1.1 -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len) + " | " + samtools_binary + " view -@ 5 -h -b -f4 -o - -U " + os.path.join(out_dir, os.path.basename(sample_file) + "_" + index_dir.split(os.sep)[-1] + "_aligned.bam")
                else:
                    cmd = snap_binary + " single " + index_dir + " -" + format + " - -o -" + format + " - -f -xf 1.1 -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len) + " | " + samtools_binary + " view -h -f4 -o - -U " + os.path.join(out_dir, os.path.basename(sample_file) + "_" + index_dir.split(os.sep)[-1] + "_aligned.sam")
            else:
                raise Exception("aligned must be 'aligned', 'unaligned', or 'all'")
            big_cmd = big_cmd + " | " + cmd
        else:
            if filt == 'aligned':
                cmd = snap_binary + " single " + index_dir + " -" + format + " - -o " + "-" + format + " " + out_file + " -f -F a -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len)
            elif filt == 'unaligned':
                cmd = snap_binary + " single " + index_dir + " -" + format + " - -o " + "-" + format + " " + out_file + " -f -F u -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len)
            elif filt == 'all':
                if format == "bam":
                    cmd = snap_binary + " single " + index_dir + " -" + format + " - -o -" + format + " - -f -xf 1.1 -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len) + " | " + samtools_binary + " view -@ 5 -h -b -f4 -o " + out_file + " -U " + os.path.join(out_dir, os.path.basename(sample_file) + "_" + index_dir.split(os.sep)[-1] + "_aligned.bam")
                else:
                    cmd = snap_binary + " single " + index_dir + " -" + format + " - -o -" + format + " - -f -xf 1.1 -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len) + " | " + samtools_binary + " view -h -f4 -o " + out_file + " -U " + os.path.join(out_dir, os.path.basename(sample_file) + "_" + index_dir.split(os.sep)[-1] + "_aligned.sam")
            else:
                raise Exception("aligned must be 'aligned', 'unaligned', or 'all'")
            big_cmd = big_cmd + " | " + cmd
        i += 1
    if len(big_cmd) >= 2616670:
        raise Exception("The typical maximum command length is 2616670, and running it with this many indicies would exceed that. Please iterate over index_dirs in chunks.")
    else:
        # exit_code = subprocess.call(big_cmd, shell=True,  stdout=FNULL, stderr=subprocess.STDOUT)
        big_cmd = "set -o pipefail; " + big_cmd
        # print(big_cmd)
        exit_code = subprocess.call(big_cmd, shell=True,  stdout=FNULL, stderr=out_message_file)
        if exit_code != 0:
            raise Exception("stream_align_single failed. Due to how snap-align prints its error messages, you will have to go digging in the file " + out_message_file_name + " to find the error. If you find an error regarding -xf, increase it slightly (say, 1.2) and try again. If you get an mmap error, you will need to decrease -xf or else try with fewer index_dirs.")
        # print(exit_code)
    # return exit_code


def sam2fastq(sam_file, out_file):
    """
    Converts a SAM file to a FASTQ file. Quite hacky
    :param sam_file: input SAM file
    :param out_file: output FASTQ file
    :return: exit code
    """
    FNULL = open(os.devnull, 'w')
    cmd = "cat " + sam_file + " | grep -v ^@ | awk '{print \"@\"$1\"\\n\"$10\"\\n+\\n\"$11}' > " + out_file  # grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}'
    exit_code = subprocess.call(cmd, shell=True,  stdout=FNULL, stderr=subprocess.STDOUT)
    FNULL.close()
    return exit_code


def top_down_align(sample_file, reference_files, index_dir, out_dir, threads=multiprocessing.cpu_count(), save_aligned=True, format="bam", large_index=True, seed_size=20, edit_distance=14, min_read_len=50, binary="snap-aligner"):
    """
    This function takes an input file and recursively aligns it to the individual reference files. Reads that are aligned are saved, and unaligned reads are passed to the next step in the alignment.
    Currently this is quite inefficient since it performs the alignment twice. Could try to use samtools to partition the reads into aligned and unaligned.
    :param sample_file: Input file on which to perform the alignment
    :param reference_files: list of references files to align agains
    :param index_dir: directory where the intermediate reference file indexes are stored
    :param out_dir: output directory for the results
    :param threads: number of threads to use
    :param save_aligned: True/False to save/not save the intermediate aligned reads
    :param format: one of 'sam' or 'bam'
    :param large_index: makes the index about 30% bigger, but faster for quick/inaccurate alignements
    :param seed_size: for initial match to begin alignment
    :param edit_distance: the maximum edit distance to allow for an alignment
    :param min_read_len: Specify the minimum read length to align, reads shorter than this (after clipping) stay unaligned.
     This should be a good bit bigger than the seed length or you might get some questionable alignments.  Default 50
    :param binary: location of the snap-aligner binary
    :return:
    """
    if format != "sam" and format != "bam":
        raise Exception("format must be one of 'sam' or 'bam'")
    sam_out_file_prev = os.path.join(out_dir, os.path.basename(sample_file) + "_unaligned_prev." + format)
    sam_out_file = os.path.join(out_dir, os.path.basename(sample_file) + "_unaligned." + format)
    for i in range(len(reference_files)):
        reference_file = reference_files[i]
        if i == 0:
            build_reference(reference_file, index_dir, large_index=large_index, seed_size=seed_size, threads=threads, binary=binary)
            if save_aligned:
                aligned_out_file = os.path.join(out_dir, os.path.basename(sample_file) + "_" + os.path.basename(reference_file) + "_" + "aligned." + format)
                align_reads(index_dir, sample_file, aligned_out_file, filt="aligned", threads=threads, edit_distance=edit_distance, min_read_len=min_read_len, binary=binary)
            align_reads(index_dir, sample_file, sam_out_file_prev, filt="unaligned", threads=threads, edit_distance=edit_distance, min_read_len=min_read_len, binary=binary)
        else:
            build_reference(reference_file, index_dir, large_index=large_index, seed_size=seed_size, threads=threads, binary=binary)
            if save_aligned:
                aligned_out_file = os.path.join(out_dir, os.path.basename(sample_file) + "_" + os.path.basename(reference_file) + "_" + "aligned." + format)
                align_reads(index_dir, sam_out_file_prev, aligned_out_file, filt="aligned", threads=threads, edit_distance=edit_distance, min_read_len=min_read_len, binary=binary)
            align_reads(index_dir, sam_out_file_prev, sam_out_file, filt="unaligned", threads=threads, edit_distance=edit_distance, min_read_len=min_read_len, binary=binary)
            shutil.move(sam_out_file, sam_out_file_prev)  # need to check if this is ok
    shutil.move(sam_out_file_prev, sam_out_file)
    return sam_out_file

def top_down_align2(sample_file, reference_files, index_dir, out_dir, threads=multiprocessing.cpu_count(), save_aligned=True, format="bam", large_index=True, seed_size=20, edit_distance=14, min_read_len=50, xf=1.2, snap_binary="snap-aligner", samtools_binary="/local/cluster/samtools/bin/samtools"):
    """
    This function takes an input file and recursively aligns it to the individual reference files. Reads that are aligned are saved, and unaligned reads are passed to the next step in the alignment.
    Currently this is quite inefficient since it performs the alignment twice. Could try to use samtools to partition the reads into aligned and unaligned.
    :param sample_file: Input file on which to perform the alignment
    :param reference_files: list of references files to align agains
    :param index_dir: directory where the intermediate reference file indexes are stored
    :param out_dir: output directory for the results
    :param threads: number of threads to use
    :param save_aligned: True/False to save/not save the intermediate aligned reads
    :param format: one of 'sam' or 'bam'
    :param large_index: makes the index about 30% bigger, but faster for quick/inaccurate alignements
    :param seed_size: for initial match to begin alignment
    :param edit_distance: the maximum edit distance to allow for an alignment
    :param min_read_len: Specify the minimum read length to align, reads shorter than this (after clipping) stay unaligned.
     This should be a good bit bigger than the seed length or you might get some questionable alignments.  Default 50
    :param binary: location of the snap-aligner binary
    :return:
    """
    FNULL = open(os.devnull, 'w')
    out_message_file_name = os.path.join(out_dir, "top_down_align_stdout.txt")
    out_message_file = open(out_message_file_name, 'w')
    if format != "sam" and format != "bam":
        raise Exception("format must be one of 'sam' or 'bam'")
    unaligned_prev = os.path.join(out_dir, os.path.basename(sample_file) + "_unaligned_prev." + format)
    unaligned = os.path.join(out_dir, os.path.basename(sample_file) + "_unaligned." + format)
    for i in range(len(reference_files)):
        reference_file = reference_files[i]
        if i == 0:
            build_reference(reference_file, index_dir, large_index=large_index, seed_size=seed_size, threads=threads, binary=snap_binary)
            aligned_out_file = os.path.join(out_dir, os.path.basename(sample_file) + "_" + os.path.basename(reference_file) + "_" + "aligned." + format)
            #cmd = "set -o pipefail; " + snap_binary + " single " + index_dir + " " + sample_file + " -o -sam - -f -xf " + str(xf) + " -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len) + " | " + samtools_binary + " view -@ " + str(threads) + " -h -b -1 -f4 -o " + unaligned_prev + " -U " + aligned_out_file
            cmd = "set -o pipefail; " + snap_binary + " single " + index_dir + " " + sample_file + " -o -sam - -xf " + str(xf) + " -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len) + " | " + samtools_binary + " view -@ " + str(threads) + " -h -b -1 -f4 -o " + unaligned_prev + " -U " + aligned_out_file
            exit_code = subprocess.check_call(cmd, shell=True,  stdout=FNULL, stderr=out_message_file)
            # align_reads(index_dir, sample_file, aligned_out_file, filt="aligned", threads=threads, edit_distance=edit_distance, min_read_len=min_read_len, binary=binary)
            # align_reads(index_dir, sample_file, sam_out_file_prev, filt="unaligned", threads=threads, edit_distance=edit_distance, min_read_len=min_read_len, binary=binary)
        else:
            build_reference(reference_file, index_dir, large_index=large_index, seed_size=seed_size, threads=threads, binary=snap_binary)
            aligned_out_file = os.path.join(out_dir, os.path.basename(sample_file) + "_" + os.path.basename(reference_file) + "_" + "aligned." + format)
            #cmd = "set -o pipefail; " + snap_binary + " single " + index_dir + " " + unaligned_prev + " -o -sam - -f -xf " + str(xf) + " -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len) + " | " + samtools_binary + " view -@ " + str(threads) + " -h -b -1 -f4 -o " + unaligned + " -U " + aligned_out_file
            #cmd = "set -o pipefail; " + samtools_binary + " view -h " + unaligned_prev + " | " + snap_binary + " single " + index_dir + " -sam - -o -sam - -f -xf " + str(xf) + " -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len) + " | " + samtools_binary + " view -@ " + str(threads) + " -h -b -1 -f4 -o " + unaligned + " -U " + aligned_out_file
            #cmd = "set -o pipefail; " + samtools_binary + " view -h " + unaligned_prev + " | " + snap_binary + " single " + index_dir + " -sam - -o -sam - -f -xf " + str(xf) + " -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len) + " | cut -f 1-18 | " + samtools_binary + " view -@ " + str(threads) + " -h -b -1 -f4 -o " + unaligned + " -U " + aligned_out_file
            cmd = "set -o pipefail; " + samtools_binary + " view -h " + unaligned_prev + " | " + snap_binary + " single " + index_dir + " -sam - -o -sam - -xf " + str(xf) + " -t " + str(threads) + " -d " + str(edit_distance) + " -mrl " + str(min_read_len) + " | cut -f 1-18 | " + samtools_binary + " view -@ " + str(threads) + " -h -b -1 -f4 -o " + unaligned + " -U " + aligned_out_file
            exit_code = subprocess.check_call(cmd, shell=True,  stdout=FNULL, stderr=out_message_file)
            # align_reads(index_dir, sample_file, aligned_out_file, filt="aligned", threads=threads, edit_distance=edit_distance, min_read_len=min_read_len, binary=binary)
            # align_reads(index_dir, sample_file, sam_out_file_prev, filt="unaligned", threads=threads, edit_distance=edit_distance, min_read_len=min_read_len, binary=binary)            shutil.move(sam_out_file, sam_out_file_prev)  # need to check if this is ok
            shutil.move(unaligned, unaligned_prev)
    shutil.move(unaligned_prev, unaligned)


def minia_assemble(out_dir, sample_file, reference_files, format='bam', top_down=False, minia_binary='/home/pi/koslickd/minia-2.0.3-Linux/bin/./minia', samtools_binary="/local/cluster/samtools/bin/samtools"):
    """
    For each one of the reference_files, convert bam to sam, sam to fastq, then assemble. Iterate while including contigs from previous runs
    :param out_dir:
    :param sample_file:
    :param reference_files:
    :param minia_binary:
    :return:
    """
    FNULL = open(os.devnull, 'w')
    out_message_file_name = os.path.join(out_dir, "minia_align_stdout.txt")
    out_message_file = open(out_message_file_name, 'w')
    unaligned = os.path.join(out_dir, os.path.basename(sample_file) + "_unaligned." + format)
    if top_down: #assemble one at a time, merging as you go
        for i in range(len(reference_files)):
            reference_file = reference_files[i]
            aligned_out_file = os.path.join(out_dir, os.path.basename(sample_file) + "_" + os.path.basename(reference_file) + "_" + "aligned." + format)
            if format == 'bam':
                sam_out = os.path.join(out_dir, "temp.sam")
                cmd = samtools_binary + " view " + aligned_out_file + " -o " + sam_out
                exit_code = subprocess.check_call(cmd, shell=True,  stdout=FNULL, stderr=out_message_file)
                sam2fastq(sam_out, os.path.join(out_dir, "temp.fastq"))
                files_to_assemble = [os.path.join(out_dir, "temp.fastq")]
            else:
                files_to_assemble = [aligned_out_file]
            if i != 0:  # Add the contigs from the previous run
                files_to_assemble.append(os.path.join(out_dir, "prev_contigs.fa"))
            fid = open(os.path.join(out_dir, "to_assemble.txt"), 'w')
            for file_name in files_to_assemble:
                fid.write(file_name + "\n")
            fid.close()
            # Run minia
            cmd = minia_binary + " -in " + os.path.join(out_dir, "to_assemble.txt") + " -kmer-size 31 -abundance-min 1 -out " + os.path.join(out_dir, "minia_out")
            exit_code = subprocess.check_call(cmd, shell=True)
            shutil.move(os.path.join(out_dir, "minia_out.contigs.fa"), os.path.join(out_dir, "prev_contigs.fa"))
        #Lastly, do the unassembled reads
        aligned_out_file = unaligned
        if format == 'bam':
            sam_out = os.path.join(out_dir, "temp.sam")
            cmd = samtools_binary + " view " + aligned_out_file + " -o " + sam_out
            exit_code = subprocess.check_call(cmd, shell=True,  stdout=FNULL, stderr=out_message_file)
            sam2fastq(sam_out, os.path.join(out_dir, "temp.fastq"))
            files_to_assemble = [os.path.join(out_dir, "temp.fastq")]
        else:
            sam2fastq(aligned_out_file, os.path.join(out_dir, "temp.fastq"))
            files_to_assemble = [os.path.join(out_dir, "temp.fastq")]
        files_to_assemble.append(os.path.join(out_dir, "prev_contigs.fa"))
        fid = open(os.path.join(out_dir, "to_assemble.txt"), 'w')
        for file_name in files_to_assemble:
            fid.write(file_name + "\n")
        fid.close()
        cmd = minia_binary + " -in " + os.path.join(out_dir, "to_assemble.txt") + " -kmer-size 31 -abundance-min 1 -out " + os.path.join(out_dir, "minia_out")
        exit_code = subprocess.check_call(cmd, shell=True)
        shutil.move(os.path.join(out_dir, "minia_out.contigs.fa"), os.path.join(out_dir, "final_assembly.fa"))
    else:  # else assemble everything individually, then join
        for i in range(len(reference_files) + 1):
            if i == 0: #Do the unaligned reads
                reference_file = unaligned
                aligned_out_file = unaligned
            else:
                reference_file = reference_files[i - 1]
                aligned_out_file = os.path.join(out_dir, os.path.basename(sample_file) + "_" + os.path.basename(reference_file) + "_" + "aligned." + format)
            if format == 'bam':
                sam_out = os.path.join(out_dir, "temp.sam")
                cmd = samtools_binary + " view " + aligned_out_file + " -o " + sam_out
                exit_code = subprocess.check_call(cmd, shell=True,  stdout=FNULL, stderr=out_message_file)
                sam2fastq(sam_out, os.path.join(out_dir, "temp.fastq"))
                file_to_assemble = os.path.join(out_dir, "temp.fastq")
            else:
                sam2fastq(aligned_out_file, os.path.join(out_dir, "temp.fastq"))
                file_to_assemble = os.path.join(out_dir, "temp.fastq")
            # Run minia
            cmd = minia_binary + " -in " + file_to_assemble + " -kmer-size 31 -abundance-min 1 -out " + os.path.join(out_dir, os.path.basename(reference_file) + "_minia_out")
            exit_code = subprocess.check_call(cmd, shell=True)
        # merge everything in one big assembly
        fid = open(os.path.join(out_dir, "to_assemble.txt"), 'w')
        for i in range(len(reference_files) + 1):
            if i == 0:
                reference_file = unaligned
            else:
                reference_file = reference_files[i - 1]
            file_name = os.path.join(out_dir, os.path.basename(reference_file) + "_minia_out.contigs.fa")
            if os.path.exists(file_name):
                fid.write(file_name + "\n")
        fid.close()
        cmd = minia_binary + " -in " + os.path.join(out_dir, "to_assemble.txt") + " -kmer-size 31 -abundance-min 1 -out " + os.path.join(out_dir, "final_assembly")
        exit_code = subprocess.check_call(cmd, shell=True)

def bam2consensus(in_bam, out_dir, gap_allow=50, samtools_binary="/local/cluster/samtools/bin/samtools"):
    out_file = os.path.join(out_dir, os.path.basename(in_bam) + ".contigs.fa")
    #out_file = '/scratch/temp/SNAP/training/test.contigs.fa'
    # run mpileup on the input bam
    cmd = samtools_binary + " sort " + in_bam + " -o " + os.path.join(out_dir, os.path.basename(in_bam)+".sorted.bam")
    exit_code = subprocess.check_call(cmd, shell=True)
    cmd = samtools_binary + " index " + os.path.join(out_dir, os.path.basename(in_bam)+".sorted.bam")
    exit_code = subprocess.check_call(cmd, shell=True)
    cmd = samtools_binary + " mpileup " + os.path.join(out_dir, os.path.basename(in_bam)+".sorted.bam") + " -o " + os.path.join(out_dir, os.path.basename(in_bam) + ".pileup")
    exit_code = subprocess.check_call(cmd, shell=True)
    os.remove(os.path.join(out_dir, os.path.basename(in_bam)+".sorted.bam"))
    os.remove(os.path.join(out_dir, os.path.basename(in_bam)+".sorted.bam.bai"))
    in_file = os.path.join(out_dir, os.path.basename(in_bam) + ".pileup")
    out_fid = open(out_file, 'w')
    fid = open(in_file, 'r')
    i = 0
    regex = re.compile(r'[actgACTG\*]')
    for line in fid.readlines():
        call = ''
        line_split = line.split()
        # get number of reads aligning to this base
        num_reads = float(line_split[3])
        # genomic position of this base
        position = int(line_split[1])
        if i == 0:
            prev_position = position - 1
            out_fid.write('>seq0\n')
            i = 1
        # If there is a small enough gap, fill with N's
        if position - prev_position <= gap_allow and position - prev_position != 1:
            for dummy in range(position - prev_position - 1):
                #out_fid.write('N')
                out_fid.write(random.choice('ACTG'))
        elif position - prev_position > gap_allow:
            out_fid.write(call)
            out_fid.write('\n')
            out_fid.write('>seq%i\n' % i)
            i += 1
        if num_reads > 0:
            chars = line_split[4]
            # get only the ACTG's
            chars_only_bases = "".join(regex.findall(chars)).upper()
            # get the most frequent
            letters_and_counts = collections.Counter(chars_only_bases).most_common()
            to_pick = []
            to_pick.append(letters_and_counts[0][0])
            for j in range(1,len(letters_and_counts)):
                if letters_and_counts[j][1] == letters_and_counts[j - 1][1]:
                    to_pick.append(letters_and_counts[j][0])
            if len(to_pick) > 1:  # If more than one maximum likelihood base, randomly pick one
                choice = random.choice(to_pick)
                if choice == '*':  # If the most common is a deletion, don't include an N
                    #call = call + ''
                    whatevs = 0
                else:
                    out_fid.write(choice)
            elif len(to_pick) == 1:  # If there's only a single maximum likelihood base, make this it
                if to_pick[0] == '*':
                    whatevs = 0
                else:
                    #call = call + to_pick[0]
                    out_fid.write(to_pick[0])
            else:  # If there is no base, make it an N
                #out_fid.write('N')
                out_fid.write(random.choice('ACTG'))
            if letters_and_counts[0][1] < num_reads / float(2):  # In case the majority of reads say there should be a deletion here
                whatevs = 0
        else:
            #out_fid.write('N')
            out_fid.write(random.choice('ACTG'))
        prev_position = position

    out_fid.write('\n')
    out_fid.close()
    fid.close()
    os.remove(os.path.join(out_dir, os.path.basename(in_bam) + ".pileup"))



##########################################################################
# Tests

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
    CE1 = CountEstimator(n=5, max_prime=9999999999971., ksize=1)
    CE2 = CountEstimator(n=5, max_prime=9999999999971., ksize=1)
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
    hash_list = {424517919, 660397082}
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
    Y = CE1.count_vector([CE1, CE2, CE3])
    assert (Y == np.array([1.,1.,0.5625])).all()
    Y2 = CE1.jaccard_vector([CE1, CE2, CE3])
    assert (Y2 == np.array([1.,1.,0.4])).all()


def form_matrices_test():
    CE1 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y')
    CE2 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y')
    CE3 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y')
    seq1 = 'tacgactgatgcatgatcgaactgatgcactcgtgatgc'
    seq2 = 'tacgactgatgcatgatcgaactgatgcactcgtgatgc'
    seq3 = 'ttgatactcaatccgcatgcatgcatgacgatgcatgatgtacgactgatgcatgatcgaactgatgcactcgtgatgczxerqwewdfhg'
    CE1.add_sequence(seq1)
    CE2.add_sequence(seq2)
    CE3.add_sequence(seq3)
    A = form_jaccard_count_matrix([CE1, CE2, CE3])
    assert (A == np.array([[1., 1., 0.80952380952380953], [1., 1., 0.80952380952380953], [0.5625, 0.5625, 1.]])).all()
    B = form_jaccard_matrix([CE1, CE2, CE3])
    assert (B == np.array([[1., 1., 0.4], [1., 1., 0.4], [0.4, 0.4, 1.]])).all()


def test_lsqnonneg():
    CE1 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y')
    CE2 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y')
    CE3 = CountEstimator(n=5, max_prime=1e10, ksize=3, save_kmers='y')
    seq1 = 'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa'
    seq2 = 'ccccccccccccccccccccccccccccccccccccccc'
    seq3 = 'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaccccccccccccccccccccccccccccccccccccccc'
    CE1.add_sequence(seq1)
    CE2.add_sequence(seq2)
    CE3.add_sequence(seq3)
    A = form_jaccard_count_matrix([CE1, CE2, CE3])
    Y = CE3.count_vector([CE1, CE2, CE3])
    x = lsqnonneg(A, Y, 0)
    # check if the L1 norm to the truth is less than machine zero
    assert abs(x - np.array([0, 0, 1.])).sum() < 1e-07
    assert (x == jaccard_count_lsqnonneg([CE1, CE2, CE3], Y, 0)).all()
    A = form_jaccard_matrix([CE1, CE2, CE3])
    Y = CE3.jaccard_vector([CE1, CE2, CE3])
    x = lsqnonneg(A, Y, 0)
    # check if the L1 norm to the truth is less than machine zero
    assert abs(x - np.array([0, 0, 1.])).sum() < 1e-07
    assert (x == jaccard_lsqnonneg([CE1, CE2, CE3], Y, 0)).all()


def test_snap():
    temp_dir = tempfile.mkdtemp()
    index_file = os.path.join(temp_dir,"in.fasta")
    align_file = os.path.join(temp_dir, "in.fastq")
    out_sam = os.path.join(temp_dir,"out.sam")
    out_fastq = os.path.join(temp_dir, "out.fastq")
    fid = open(index_file, "w")
    fid.write(">seq1\n")
    fid.write("GGATTGGTGTATTCACGCTAGAATTCTTGTTAATCATATTATAACACTGGTTAATAGAGGAATGCAAAAAGATGC\n")
    fid.close()
    fid = open(align_file, "w")
    fid.write("@SRR172902.1213325\n")
    fid.write("GGATTGGTGTATTCACGCTAGAATTCTTGTTAATCATATTATAACACTGGTTAATAGAGGAATGCAAAAAGATGC\n")
    fid.write("+\n")
    fid.write("BB@BABB<B>BBBBBBBBABB@@AAAA@A<@@@@A@A@@?@A?A@?@AB@;@@>@@A?>@>@<?@?=9==>?<<@\n")
    fid.close()
    res = build_reference(index_file, temp_dir, large_index=True, seed_size=20, threads=1, binary="snap-aligner")
    assert res == 0
    res = align_reads(temp_dir, align_file, out_sam, filt='aligned', threads=1, edit_distance=14, min_read_len=50, binary="snap-aligner")
    assert res == 0
    sam2fastq(out_sam, out_fastq)
    assert filecmp.cmp(align_file, out_fastq)
    shutil.rmtree(temp_dir)


def test_suite():
    """
    Runs all the test functions
    :return:
    """
    from sys import platform as _platform
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
    form_matrices_test()
    test_lsqnonneg()
    if _platform == "linux" or _platform == "linux2":  # I don't have snap installed on my mac, so only do the test on the server
        test_snap()

