import MinHash as MH

reload(MH)
fid = open('/Users/dkoslicki/Dropbox/Repositories/MinHash/data/test_files.txt', 'r')
file_names = fid.readlines()
fid.close()
file_names = [name.strip() for name in file_names]
CE = MH.CountEstimator(n=500, ksize=11, input_file_name=file_names[0], save_kmers='y')
CE2 = MH.CountEstimator(n=500, ksize=11, input_file_name=file_names[1], save_kmers='y')
CE2.jaccard_count(CE)
CE2.jaccard(CE)
CE.jaccard(CE2)

CS = MH.CompositionSketch(n=5000, ksize=11, prefixsize=1, input_file_name=file_names[0])
CS2 = MH.CompositionSketch(n=5000, ksize=11, prefixsize=1, input_file_name=file_names[1])
CS.jaccard_count(CS2)
CS.jaccard(CS2)
CS2.jaccard(CS)


i = 0
for record in screed.open(file_names[0]):
    for kmer in MH.kmers(record.sequence,10):
        if kmer == 'TGGAATTCCA':
            i += 1

print(i)

CE = MH.CountEstimator(n=500, ksize=20, input_file_name=file_names[0], save_kmers='y')

#Compute a bunch of them and save to a single, big HDF5 file
CEs = MH.compute_multiple(n=500,ksize=11,input_files_list=file_names,save_kmers='y')
MH.export_multiple_to_single_hdf5(CEs, 'test_big.h5')
#Load them back in
CEs=MH.import_multiple_from_single_hdf5('test_big.h5')
#load just a few in
CEs=MH.import_multiple_from_single_hdf5('test_big.h5', file_names[0:2])


# Let's look at forming the Y vector
CEs = MH.import_multiple_hdf5(out_file_names)
MCE = MH.import_single_hdf5('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N500_k31_inComparison.h5')
Y = np.zeros(len(CEs))
for i in range(len(CEs)):
    Y[i] = CEs[i].jaccard_count(MCE)[1]



################################
# Old, no chunking
def form_common_kmer_matrix(CEs):
    """
    Forms the common kmer matrix for input list of count estimators
    :param CEs: a list of Count Estimators
    :return: a numpy matrix A where A_{i,j} \approx \sum_{w\in SW_k(g_i) \cap SW_k{g_j}} \frac{occ_w(g_j)}{|g_j| - k + 1}
    """
    # I could decreae the memory usage by not creating the whole list of input_args. Use an exterior chunk_size to create
    # a smaller list of CEs, consume those, and then repeat until finished.
    A = np.zeros((len(CEs), len(CEs)), dtype=np.float64)
    input_args = collections.deque()
    indicies = []
    for i in xrange(len(CEs)):
        for j in xrange(len(CEs)):
            input_args.append((CEs[i], CEs[j]))
            indicies.append((i, j))

    #input_args = ((CEs[i], CEs[j]) for i in xrange(len(CEs)) for j in xrange(len(CEs)))
    pool = Pool(processes=multiprocessing.cpu_count())
    res = pool.imap(form_common_kmer_matrix_helper, input_args, chunksize=np.floor(len(indicies)/float(multiprocessing.cpu_count())))  # chunk into fewest pieces possible
    # pool.close()
    # pool.join()
    # pool.terminate()
    #for i in xrange(len(indicies)):
    for i, val in enumerate(res):
        A[indicies[i][0], indicies[i][1]] = val[0] #res[i][0]  # Replace i with i+last_index where last_index was the number of times the xranges executed before going into the pool
        A[indicies[i][1], indicies[i][0]] = val[1] #res[i][1]

    pool.terminate()
    return A


##############################
# Try exterior chunking, still takes to long to form the input_args
def form_common_kmer_matrix(all_CEs):
    """
    Forms the common kmer matrix for input list of count estimators
    :param all_CEs: a list of Count Estimators
    :return: a numpy matrix A where A_{i,j} \approx \sum_{w\in SW_k(g_i) \cap SW_k{g_j}} \frac{occ_w(g_j)}{|g_j| - k + 1}
    """
    A = np.zeros((len(all_CEs), len(all_CEs)), dtype=np.float64)
    chunk_size = 70000
    # Can precompute all the indicies
    indicies = []
    for i in xrange(len(all_CEs)):
        for j in xrange(len(all_CEs)):
            indicies.append((i, j))
    for sub_indicies in chunks(indicies, chunk_size):
        input_args = ((all_CEs[sub_indicies[i][0]], all_CEs[sub_indicies[i][1]]) for i in xrange(len(sub_indicies)))
        pool = Pool(processes=multiprocessing.cpu_count())
        res = pool.imap(form_common_kmer_matrix_helper, input_args, chunksize=np.floor(len(indicies)/float(multiprocessing.cpu_count())))
        # pool.close()
        # pool.join()
        # pool.terminate()
        for (i, j), val in zip(sub_indicies, res):
            A[i, j] = val[0] #res[i][0]  # Replace i with i+last_index where last_index was the number of times the xranges executed before going into the pool
            A[j, i] = val[1] #res[i][1]
            print((i,j))
            print(val)

        pool.terminate()
    return A


############################
# Try with shared array
def temp_func(ij,shared_mins=shared_mins, shared_counts=shared_counts, p=p):
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
    for (count1, count2) in yield_count_overlaps(mins1, mins2, counts1, counts2):
        common1 += count1  # The unpopulated hashes have count 0, so we don't have to worry about that here
        common2 += count2
    return (common2 / float(sum(counts2)), common1 / float(sum(counts1)))

import multiprocessing
import ctypes
import numpy as np

def form_common_kmer_matrix(all_CEs):
    A = np.zeros((len(all_CEs), len(all_CEs)), dtype=np.float64)
    # Can precompute all the indicies
    indicies = []
    for i in xrange(len(all_CEs)):
        for j in xrange(len(all_CEs)):
            indicies.append((i, j))

    shared_mins_base = multiprocessing.Array(ctypes.c_longlong, len(all_CEs)*len(all_CEs[0]._mins))
    shared_mins = np.ctypeslib.as_array(shared_mins_base.get_obj())
    shared_mins = shared_mins.reshape(len(all_CEs), len(all_CEs[0]._mins))
    shared_counts_base = multiprocessing.Array(ctypes.c_double, len(all_CEs)*len(all_CEs[0]._counts))
    shared_counts = np.ctypeslib.as_array(shared_counts_base.get_obj())
    shared_counts = shared_counts.reshape(len(all_CEs), len(all_CEs[0]._counts))
    p = all_CEs[0].p
    for i in range(len(all_CEs)):
        shared_mins[i] = all_CEs[i]._mins
        shared_counts[i] = all_CEs[i]._counts

    pool = multiprocessing.Pool(processes=4)
    res = pool.map(temp_func, indicies)
    for (i, j), val in zip(indicies, res):
        A[i, j] = val[0] #res[i][0]  # Replace i with i+last_index where last_index was the number of times the xranges executed before going into the pool
        A[j, i] = val[1] #res[i][1]
        print((i,j))
        print(val)

    pool.terminate()
    return A