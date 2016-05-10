import sys
sys.path.append('/nfs1/Koslicki_Lab/koslickd/Repositories/MinHashMetagenomics/src/')
import os, timeit, h5py
import MinHash as MH
import numpy as np

fid = open('/nfs1/Koslicki_Lab/koslickd/MinHash/Data/FileNames.txt', 'r')
file_names = fid.readlines()
fid.close()
file_names = [name.strip() for name in file_names]

###############################
# Compute the hashes for all the training genomes
n = 500
CEs = MH.compute_multiple(n=n, max_prime=9999999999971., ksize=31, input_files_list=file_names, save_kmers='y', num_threads=48)
# Export
MH.export_multiple_hdf5(CEs, '/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N'+str(n)+'k31/')
# Save the hashes in the training genomes
hash_list = set()
for CE in CEs:
    hash_list.update(CE._mins)

fid = h5py.File('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N'+str(n)+'k31_mins.h5', 'w')
fid.create_dataset("hash_list", data=list(hash_list))
fid.close()
# If I need to read it back in
#fid = h5py.File('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N500k31_mins.h5','r')
#hash_list = set(fid["hash_list"][:])

n = 5000
CEs = MH.compute_multiple(n=n, max_prime=9999999999971., ksize=31, input_files_list=file_names, save_kmers='y', num_threads=48)
# Export
MH.export_multiple_hdf5(CEs, '/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N'+str(n)+'k31/')
# Save the hashes in the training genomes
hash_list = set()
for CE in CEs:
    hash_list.update(CE._mins)

fid = h5py.File('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N'+str(n)+'k31_mins.h5', 'w')
fid.create_dataset("hash_list", data=list(hash_list))
fid.close()
# If I need to read it back in
#fid = h5py.File('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N500k31_mins.h5','r')
#hash_list = set(fid["hash_list"][:])

n = 50000
hash_list = set()
chunk_size = 500
for i in range(0, len(file_names), chunk_size):
    temp_file_names = file_names[i:i+chunk_size]
    CEs = MH.compute_multiple(n=n, max_prime=9999999999971., ksize=31, input_files_list=temp_file_names, save_kmers='y', num_threads=48)
    # Export
    MH.export_multiple_hdf5(CEs, '/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N'+str(n)+'k31/')
    # Save the hashes in the training genomes
    for CE in CEs:
        hash_list.update(CE._mins)

fid = h5py.File('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N'+str(n)+'k31_mins.h5', 'w')
fid.create_dataset("hash_list", data=list(hash_list))
fid.close()
# If I need to read it back in
#fid = h5py.File('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N500k31_mins.h5','r')
#hash_list = set(fid["hash_list"][:])

####################################
# Form a CE for a metagenome
n = 500
fid = h5py.File('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N'+str(n)+'k31_mins.h5', 'r')
hash_list = set(fid["hash_list"][...])
fid.close()
CE = MH.CountEstimator(n=n, max_prime=9999999999971., ksize=31, input_file_name='/nfs1/Koslicki_Lab/koslickd/MinHash/Data/SRR172902.fastq', save_kmers='y')
CE.export('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N'+str(n)+'_k31_all.h5')
CE2 = MH.CountEstimator(n=n, max_prime=9999999999971., ksize=31, input_file_name='/nfs1/Koslicki_Lab/koslickd/MinHash/Data/SRR172902.fastq', save_kmers='y', hash_list=hash_list)
CE2.export('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N'+str(n)+'_k31_inComparison.h5')
n = 5000
fid = h5py.File('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N'+str(n)+'k31_mins.h5', 'r')
hash_list = set(fid["hash_list"][...])
fid.close()
MCE2 = MH.CountEstimator(n=n, max_prime=9999999999971., ksize=31, input_file_name='/nfs1/Koslicki_Lab/koslickd/MinHash/Data/SRR172902.fastq', save_kmers='y', hash_list=hash_list)
MCE2.export('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N'+str(n)+'_k31_inComparison.h5')
MCE = MH.CountEstimator(n=n, max_prime=9999999999971., ksize=31, input_file_name='/nfs1/Koslicki_Lab/koslickd/MinHash/Data/SRR172902.fastq', save_kmers='y')
MCE.export('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N'+str(n)+'_k31_all.h5')
n = 50000
if os.path.isfile('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N'+str(n)+'k31_mins.h5'):
    fid = h5py.File('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N'+str(n)+'k31_mins.h5', 'r')
    hash_list = set(fid["hash_list"][...])
    fid.close()
    MCE2 = MH.CountEstimator(n=n, max_prime=9999999999971., ksize=31, input_file_name='/nfs1/Koslicki_Lab/koslickd/MinHash/Data/SRR172902.fastq', save_kmers='y', hash_list=hash_list)
    MCE2.export('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N'+str(n)+'_k31_inComparison.h5')
    MCE = MH.CountEstimator(n=n, max_prime=9999999999971., ksize=31, input_file_name='/nfs1/Koslicki_Lab/koslickd/MinHash/Data/SRR172902.fastq', save_kmers='y')
    MCE.export('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N'+str(n)+'_k31_all.h5')

exit()

###################################
# Make the Y vectors
# Read in all the saved hashes
import sys
sys.path.append('/nfs1/Koslicki_Lab/koslickd/Repositories/MinHashMetagenomics/src/')
import os, timeit, h5py
import MinHash as MH
import numpy as np
import logging
fid = open('/nfs1/Koslicki_Lab/koslickd/MinHash/Data/FileNames.txt', 'r')
file_names = fid.readlines()
fid.close()
file_names = [name.strip() for name in file_names]
training_n = 5000
out_file_names = ["/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N"+str(training_n)+"k31/" + os.path.basename(item) + ".CE.h5" for item in file_names]
CEs = MH.import_multiple_hdf5(out_file_names)
# Set which CE for the test metagenome to read in
n = 50000
MCE_in_comparison = MH.import_single_hdf5('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N'+str(n)+'_k31_inComparison.h5')
MCE_all = MH.import_single_hdf5('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N'+str(n)+'_k31_all.h5')
Y_count_in_comparison = MCE_in_comparison.count_vector(CEs)
Y_jaccard_in_comparison = MCE_in_comparison.jaccard_vector(CEs)
Y_count_all = MCE_all.count_vector(CEs)
Y_jaccard_all = MCE_all.jaccard_vector(CEs)
# Read in the taxonomy and see which are the largest entries of which Y vectors
fid = open('/nfs1/Koslicki_Lab/koslickd/MinHash/Data/Taxonomy.txt', 'r')
taxonomy = fid.readlines()
fid.close()
taxonomy = [item.strip() for item in taxonomy]
taxonomy_names = [item.split('\t')[0] for item in taxonomy]

for test_Y in [Y_count_in_comparison, Y_jaccard_in_comparison, Y_count_all, Y_jaccard_all]:
    i = 0
    print("Y_vales")
    for pair in sorted(enumerate(test_Y), key=lambda x: x[1])[::-1]:
        index = pair[0]
        value = pair[1]
        print("\t name: %s abundance: %f" %(taxonomy_names[index], value))
        i += 1
        if i > 20:
            break

# n=500, Y_jaccard_in_comparison did quite well
# n=50,000 is definitely the best so far.
# Let's see if increasing the n for the training genomes helps as well
# Yes indeed it does. trainging n=5000 and sample n=50,000 gives great results, especially for the jaccard_count
# though the curious thing is that jaccard works better than jaccard_count when the n's are smaller


#################################
# Make the common Kmer matrix

import sys
sys.path.append('/nfs1/Koslicki_Lab/koslickd/Repositories/MinHashMetagenomics/src/')
import os, timeit, h5py
import MinHash as MH
import numpy as np
import logging
fid = open('/nfs1/Koslicki_Lab/koslickd/MinHash/Data/FileNames.txt', 'r')
file_names = fid.readlines()
fid.close()
file_names = [name.strip() for name in file_names]
training_n = 50000
out_file_names = ["/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N"+str(training_n)+"k31/" + os.path.basename(item) + ".CE.h5" for item in file_names]
CEs = MH.import_multiple_hdf5(out_file_names)

A = MH.form_jaccard_count_matrix(CEs)  #NOTE!!! I only need to form this for the indicies where Y[i] > 0


################################
# Test the lsqnonneg stuff
# Read in all the saved hashes
import sys, os
sys.path.append('/nfs1/Koslicki_Lab/koslickd/Repositories/MinHashMetagenomics/src/')
import MinHash as MH
import numpy as np
fid = open('/nfs1/Koslicki_Lab/koslickd/MinHash/Data/FileNames.txt', 'r')
file_names = fid.readlines()
fid.close()
file_names = [name.strip() for name in file_names]
training_n = 50000
out_file_names = ["/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N"+str(training_n)+"k31/" + os.path.basename(item) + ".CE.h5" for item in file_names]
CEs = MH.import_multiple_hdf5(out_file_names)
# Set which CE for the test metagenome to read in
n = 50000
MCE_in_comparison = MH.import_single_hdf5('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N'+str(n)+'_k31_inComparison.h5')
MCE_all = MH.import_single_hdf5('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N'+str(n)+'_k31_all.h5')
Y_count_in_comparison = MCE_in_comparison.count_vector(CEs)
Y_jaccard_in_comparison = MCE_in_comparison.jaccard_vector(CEs)
Y_count_all = MCE_all.count_vector(CEs)
Y_jaccard_all = MCE_all.jaccard_vector(CEs)
# Read in the taxonomy and see which are the largest entries of which Y vectors
fid = open('/nfs1/Koslicki_Lab/koslickd/MinHash/Data/Taxonomy.txt', 'r')
taxonomy = fid.readlines()
fid.close()
taxonomy = [item.strip() for item in taxonomy]
taxonomy_names = [item.split('\t')[0] for item in taxonomy]

reconstruction = MH.jaccard_count_lsqnonneg(CEs, Y_count_in_comparison, .0001)
i = 0
print("Reconstruction Values")
for pair in sorted(enumerate(reconstruction), key=lambda x: x[1])[::-1]:
    index = pair[0]
    value = pair[1]
    print("\t name: %s abundance: %f" %(taxonomy_names[index], value))
    i += 1
    if i > 20:
        break




for test_Y in [Y_count_in_comparison, Y_jaccard_in_comparison, Y_count_all, Y_jaccard_all]:
    i = 0
    print("Y_vales")
    for pair in sorted(enumerate(test_Y), key=lambda x: x[1])[::-1]:
        index = pair[0]
        value = pair[1]
        print("\t name: %s abundance: %f" %(taxonomy_names[index], value))
        i += 1
        if i > 20:
            break




#########################
# Check the stats of the hashed kmers
Acounts=list()
Ccounts=list()
Tcounts=list()
Gcounts=list()
for kmer in sorted(MCE._kmers):
    Acount=0
    Ccount=0
    Tcount=0
    Gcount=0
    for letter in kmer:
        if letter == 'A':
            Acount+=1
        elif letter == 'C':
            Ccount+=1
        elif letter =='T':
            Tcount+=1
        elif letter =='G':
            Gcount+=1
    Acount = Acount/float(len(kmer))
    Ccount = Ccount/float(len(kmer))
    Tcount = Tcount/float(len(kmer))
    Gcount = Gcount/float(len(kmer))
    Acounts.append(Acount)
    Ccounts.append(Ccount)
    Tcounts.append(Tcount)
    Gcounts.append(Gcount)
    print("%s\t%f\t%f\t%f\t%f" % (kmer,Acount,Ccount,Tcount,Gcount))

np.mean(Acounts)
np.mean(Ccounts)
np.mean(Tcounts)
np.mean(Gcounts)


######################
# Let's test where the longest time is spent for building the hash from a file
import sys
sys.path.append('/nfs1/Koslicki_Lab/koslickd/Repositories/MinHashMetagenomics/src/')
import os, timeit, h5py
import MinHash as MH
import numpy as np
file_name = '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000467715.fna'
t0 = timeit.default_timer()
CE = MH.CountEstimator(n=50000, max_prime=9999999999971., ksize=31, input_file_name=file_name, save_kmers='n')
t1 = timeit.default_timer()
print("Parse time: %f" % (t1-t0))


temp = [123513423523463457]*50000
t0 = timeit.default_timer()
for i in xrange(1000):
    for j in xrange(25000):
        if j == 25000:
            temp.insert(25000,1)
            temp.pop()

t1 = timeit.default_timer()
print(t1-t0)



###################
# Timing of the default methods
import timeit
import blist
import random
import numpy as np
len_list = 50000
upper = 1e15
mins = [upper]*len_list
num_its = int(1e5)
h_list = [random.randint(0, upper) for it in xrange(num_its)]
num_times = 10
times = np.zeros(num_times)

for outer_iter in xrange(num_times):
    mins = [upper]*len_list
    h_list = [random.randint(0, upper) for it in xrange(num_its)]
    t0 = timeit.default_timer()
    for h in h_list:
        if h < mins[-1]:
            i = 0
            for v in mins:  # Let's get rid of the enumerate, replace with xrange and i += 1 stuff. O.W. try blist
                if h < v:
                    mins.insert(i, h)
                    dummy = mins.pop()
                    break
                elif h == v:
                    break
                i += 1
    t1 = timeit.default_timer()
    times[outer_iter] = (t1-t0)

print("Time: %f" % np.mean(times))  # (default) enumerate, standard list. (len_list=5000, 10.3896) (len_list=50000, 343.972694)
print("Time: %f" % np.mean(times))  # in mins i+=1, standard list, (len_list=5000, 11.822465) (len_list=50000, )


# Blist
import timeit
from blist import *
import random
import numpy as np
import bisect
len_list = 50000
upper = 1e15
mins = blist([upper]*len_list)
num_its = int(1e5)
h_list = [random.randint(0, upper) for it in xrange(num_its)]
num_times = 10
times = np.zeros(num_times)

for outer_iter in xrange(num_times):
    mins = [upper]*len_list
    h_list = [random.randint(0, upper) for it in xrange(num_its)]
    t0 = timeit.default_timer()
    for h in h_list:
        if h < mins[-1]:
            i = bisect.bisect_left(mins, h)
            if mins[i] == h:  # this means that h is already in mins. This is where I would increment counts
                pass
            else:
                mins.insert(i, h)
                dummy = mins.pop()
    t1 = timeit.default_timer()
    times[outer_iter] = (t1-t0)

print("Time: %f" % np.mean(times))

