import sys
sys.path.append('/nfs1/Koslicki_Lab/koslickd/Repositories/MinHashMetagenomics/src/')
import os, timeit, h5py
import MinHash as MH
import numpy as np


fid = open('/nfs1/Koslicki_Lab/koslickd/MinHash/Data/FileNames.txt', 'r')
file_names = fid.readlines()
fid.close()
file_names = [name.strip() for name in file_names]
out_file_names = ["/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N500k31/" + os.path.basename(item) + ".CE.h5" for item in file_names]

###############################
# Compute the hashes for all the training genomes
n = 500
CEs = MH.compute_multiple(n=n, max_prime=1e10, ksize=31, input_files_list=file_names, save_kmers='y', num_threads=48)
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
CEs = MH.compute_multiple(n=n, max_prime=1e10, ksize=31, input_files_list=file_names, save_kmers='y', num_threads=48)
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
CEs = MH.compute_multiple(n=n, max_prime=1e10, ksize=31, input_files_list=file_names, save_kmers='y', num_threads=48)
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

####################################
# Form a CE for a metagenome
# n=500
CE = MH.CountEstimator(n=500, max_prime=1e10, ksize=31, input_file_name='/nfs1/Koslicki_Lab/koslickd/MinHash/Data/SRR172902.fastq', save_kmers='y')
CE.export('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N500_k31_all.h5')
CE2 = MH.CountEstimator(n=500, max_prime=1e10, ksize=31, input_file_name='/nfs1/Koslicki_Lab/koslickd/MinHash/Data/SRR172902.fastq', save_kmers='y', hash_list=hash_list)
CE2.export('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N500_k31_inComparison.h5')
# n=5000
MCE2 = MH.CountEstimator(n=5000, max_prime=1e10, ksize=31, input_file_name='/nfs1/Koslicki_Lab/koslickd/MinHash/Data/SRR172902.fastq', save_kmers='y', hash_list=hash_list)
MCE2.export('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N5000_k31_inComparison.h5')
MCE = MH.CountEstimator(n=5000, max_prime=1e10, ksize=31, input_file_name='/nfs1/Koslicki_Lab/koslickd/MinHash/Data/SRR172902.fastq', save_kmers='y')
MCE.export('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N5000_k31_all.h5')
# n=50,000
MCE2 = MH.CountEstimator(n=5000, max_prime=1e10, ksize=31, input_file_name='/nfs1/Koslicki_Lab/koslickd/MinHash/Data/SRR172902.fastq', save_kmers='y', hash_list=hash_list)
MCE2.export('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N5000_k31_inComparison.h5')
MCE = MH.CountEstimator(n=5000, max_prime=1e10, ksize=31, input_file_name='/nfs1/Koslicki_Lab/koslickd/MinHash/Data/SRR172902.fastq', save_kmers='y')
MCE.export('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N5000_k31_all.h5')

###################################
# Make the Y vectors
# Read in all the saved hashes
import sys
sys.path.append('/nfs1/Koslicki_Lab/koslickd/Repositories/MinHashMetagenomics/src/')
import os, timeit, h5py
import MinHash as MH
import numpy as np
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
Y_count_in_comparison = MCE_in_comparison.jaccard_count_vector(CEs)
Y_jaccard_in_comparison = MCE_in_comparison.jaccard_vector(CEs)
Y_count_all = MCE_all.jaccard_count_vector(CEs)
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



