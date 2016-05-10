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
