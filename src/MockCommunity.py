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
    print("Y_values")
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
# NOTE: It's probably best to precompute the A matrices, then use plain MH.lsqnonneg()
# Read in all the saved hashes
import sys, os
sys.path.append('/nfs1/Koslicki_Lab/koslickd/Repositories/MinHashMetagenomics/src/')
import MinHash as MH
import numpy as np
fid = open('/nfs1/Koslicki_Lab/koslickd/MinHash/Data/FileNames.txt', 'r')
file_names = fid.readlines()
fid.close()
file_names = [name.strip() for name in file_names]
training_n = 5000
eps = .001
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

vectors = [Y_count_in_comparison, Y_jaccard_in_comparison, Y_count_all, Y_jaccard_all]
for it in range(len(vectors)):
    test_Y = vectors[it]
    if it == 0:
        print("jaccard count in comparison")
    if it == 1:
        print("jaccard in comparison")
    if it == 2:
        print("jaccard count in all")
    if it ==3:
        print("jaccard in all")
    reconstruction = MH.jaccard_count_lsqnonneg(CEs, test_Y, eps)
    i = 0
    print("Reconstruction Values")
    for pair in sorted(enumerate(reconstruction), key=lambda x: x[1])[::-1]:
        index = pair[0]
        value = pair[1]
        print("\t name: %s abundance: %f" %(taxonomy_names[index], value))
        i += 1
        if i > 25:
            break


# Plain ol' Jaccard is looking pretty good when combined with the nnlsq!


#Let's now work on getting bowtie2 integrated, so we can do a top-down approach
Y_count_in_comparison.index(Y_count_in_comparison)
list(Y_count_in_comparison).index(max(Y_count_in_comparison))
taxonomy[1392]
###########
# bowtie2
#May first discard the reads that don't have any matched k-mers for k = argument of -L
#/local/cluster/bin/bowtie2-build -f -o 1 -t 15 /nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/G000008565.fna /nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/ref

#/local/cluster/bowtie2-2.2.3/bowtie2 --threads 48 -k 1 -R 1 --local -D 1 -L 25 -q -x ref -U /nfs1/Koslicki_Lab/koslickd/MinHash/Data/SRR172902.fastq --un unaligned.fastq --al aligned.fastq > /dev/null  # 37.67% overall alignment rate, time: 15989.11user 5797.74system 9:11.10elapsed
#/local/cluster/bowtie2-2.2.3/bowtie2 --threads 48 -k 1 -R 1 --local -D 1 -L 20 -q -x ref -U /nfs1/Koslicki_Lab/koslickd/MinHash/Data/SRR172902.fastq --un unaligned.fastq --al aligned.fastq > /dev/null  # 38.08% overall alignment rate, time: 16095.37user 5941.72system 9:17.17elapsed

###########
# CASHX
# Perhaps try CASHX as a first pass (exact matches, or up to 2bp mismatch), unfortunately requires fasta instead of fastq
#/local/cluster/CASHX_2.1/bin/./cashx_formatDB -i G000008565.fna -l 12
#/local/cluster/CASHX_2.1/bin/./cashx_formatDB

################
# SNAP
# Or try SNAP. Works only on math1
#~/./snap-aligner index /nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/G000008565.fna /nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/SNAP
#~/./snap-aligner single /nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/SNAP /nfs1/Koslicki_Lab/koslickd/MinHash/Data/SRR172902.fastq -o /nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/SNAP/out.sam -F u -t 48 -d 20
# Total Reads    Aligned, MAPQ >= 10    Aligned, MAPQ < 10     Unaligned              Too Short/Too Many Ns  %Pairs	Reads/s   Time in Aligner (s)
# 6,562,065      1,920,800 (29.27%)     69,155 (1.05%)         3,786,810 (57.71%)     785,300 (11.97%)       0.00%%	546,793   12

#/home/pi/koslickd/./snap-aligner index /nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/G000008565.fna /nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/SNAP
#/home/pi/koslickd/./snap-aligner single /nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/SNAP /nfs1/Koslicki_Lab/koslickd/MinHash/Data/SRR172902.fastq -o /nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/SNAP/out.sam -F a -t 48 -d 20
#cat out.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > aligned.fastq
#/home/pi/koslickd/./snap-aligner single /nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/SNAP /nfs1/Koslicki_Lab/koslickd/MinHash/Data/SRR172902.fastq -o /nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/SNAP/out.sam -F u -t 48 -d 20
#cat out.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > unaligned.fastq

# Build the index
import subprocess
import os
FNULL = open(os.devnull, 'w')
ind = 1392
index_dir = "/nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/SNAP"
out_dir = "/nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/"
sam_out_file = os.path.join(out_dir, "out.sam")
sample_file = "/nfs1/Koslicki_Lab/koslickd/MinHash/Data/SRR172902.fastq"
# reference_file = CEs[ind].input_file_name # Until I get snap on math0, or scipy on math1
reference_file = "/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000008565.fna"
fastq_out_basename = os.path.join(out_dir, os.path.basename(sample_file) + "_" + os.path.basename(reference_file))



MH.build_reference(reference_file, index_dir)
MH.align_reads(index_dir, sample_file, sam_out_file, filt='aligned')
MH.sam2fastq(sam_out_file, fastq_out_basename + "_aligned.fastq")
MH.align_reads(index_dir, sample_file, sam_out_file, filt='unaligned')
MH.sam2fastq(sam_out_file, fastq_out_basename + "_unaligned.fastq")
os.remove(sam_out_file)

# More efficient thing to do is figure out how to get the output SAM to split over aligned and unaligned, rather than reading it twice...
# Since snap-aligner also accepts streaming inputs, if I didn't care to save the aligned reads, I could just chain a bunch of calls to it together with different references!

###################################
# Test packaged top down approach
import sys
sys.path.append('/nfs1/Koslicki_Lab/koslickd/Repositories/MinHashMetagenomics/src/')
import os
import MinHash as MH

fid = open('/nfs1/Koslicki_Lab/koslickd/MinHash/Data/FileNames.txt', 'r')
file_names = fid.readlines()
fid.close()
file_names = [name.strip() for name in file_names]
training_n = 5000
out_file_names = ["/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N"+str(training_n)+"k31/" + os.path.basename(item) + ".CE.h5" for item in file_names]
CEs = MH.import_multiple_hdf5(out_file_names)
n = 50000
MCE_in_comparison = MH.import_single_hdf5('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N'+str(n)+'_k31_inComparison.h5')
MCE_all = MH.import_single_hdf5('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N'+str(n)+'_k31_all.h5')
Y_count_in_comparison = MCE_in_comparison.count_vector(CEs)

reference_files = []
reconstruction = MH.jaccard_count_lsqnonneg(CEs, Y_count_in_comparison, .001)
N = 40
indicies = reconstruction.argsort()[-N:][::-1]
for index in indicies:
    reference_files.append(CEs[index].input_file_name)

#index_dir = "/nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/SNAP"
index_dir = "/scratch/temp/SNAP/"
out_dir = "/nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/"
sample_file = "/nfs1/Koslicki_Lab/koslickd/MinHash/Data/SRR172902.fastq"

#reference_files = ['/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000008565.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000015425.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000012825.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000006885.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000016965.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000008525.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000012905.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000387985.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000008805.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000154225.fna']
#reference_files = ['/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000008565.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000015425.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000012825.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000006885.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000016965.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000008525.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000012905.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000387985.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000008805.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000154225.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000009005.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000007465.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000376205.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000158415.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000632785.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000016525.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000007785.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000210115.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000008285.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000248215.fna']
#reference_files = ['/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000008565.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000015425.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000012825.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000006885.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000016965.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000008525.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000012905.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000387985.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000008805.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000154225.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000009005.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000007465.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000376205.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000158415.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000632785.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000016525.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000007785.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000210115.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000008285.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000248215.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000286335.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000478485.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000012025.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000338275.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000612145.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000144085.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000003645.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000160215.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000014625.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000183905.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000025085.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000299175.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000175635.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000185885.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000309525.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000260275.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000006925.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000252445.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000696815.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000060285.fna']

out_sam = MH.top_down_align(sample_file, reference_files[0:10], index_dir, out_dir)
pre, ext = os.path.splitext(out_sam)
out_fastq = pre + ".fastq"
MH.sam2fastq(out_sam, out_fastq)


##############################
# Try the streaming approach
import sys
sys.path.append('/nfs1/Koslicki_Lab/koslickd/Repositories/MinHashMetagenomics/src/')
import os
import MinHash as MH
import shutil
import timeit
outT0 = timeit.default_timer()
fid = open('/nfs1/Koslicki_Lab/koslickd/MinHash/Data/FileNames.txt', 'r')
file_names = fid.readlines()
fid.close()
file_names = [name.strip() for name in file_names]
training_n = 50000
out_file_names = ["/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N"+str(training_n)+"k31/" + os.path.basename(item) + ".CE.h5" for item in file_names]
CEs = MH.import_multiple_hdf5(out_file_names)
n = 50000
MCE_in_comparison = MH.import_single_hdf5('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N'+str(n)+'_k31_inComparison.h5')
MCE_all = MH.import_single_hdf5('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SRR172902.fastq.CE_N'+str(n)+'_k31_all.h5')
Y_count_in_comparison = MCE_in_comparison.count_vector(CEs)
eps = .001
reconstruction = MH.jaccard_count_lsqnonneg(CEs, Y_count_in_comparison, eps)
#N = 40
# indicies = reconstruction.argsort()[-N:][::-1]
reconstruction = reconstruction/float(sum(reconstruction))
indicies = reconstruction.argsort()[::-1]
reference_files = []
for index in indicies:
    reference_files.append(CEs[index].input_file_name)
    if reconstruction[index] < eps:  # Add at leat one more index below the threshold (to make sure I don't have an empty list of indexes)
        break

print("Number of indicies to build: %d" % len(reference_files))

index_dir = "/scratch/temp/SNAP/"
out_dir = "/nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/"
sample_file = "/nfs1/Koslicki_Lab/koslickd/MinHash/Data/SRR172902.fastq"

index_dirs = MH.build_references(reference_files, index_dir)
out_sam = "/nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/out.sam"
t0 = timeit.default_timer()
MH.stream_aligned_save_unaligned(index_dirs, sample_file, out_sam)
t1 = timeit.default_timer()
print("Alignment time: %f" % (t1-t0))
pre, ext = os.path.splitext(out_sam)
out_fastq = pre + ".fastq"
MH.sam2fastq(out_sam, out_fastq)
outT1 = timeit.default_timer()
print("Total time: %f" % (outT1-outT0))

# Clean up
os.remove(out_sam)
for index_dir in index_dirs:
    shutil.rmtree(index_dir)

# Total running time, all eps = .001
# training_n 500
# total time: 172.598919
# 227M fastq

# training_n 5000
# total time: 328.972598
# 182M fastq

# training_n 50000
# total time: 624.505597
# 182M fastq




# BAM vs SAM
# training_n 5000
# 230 for SAM
# 209.065743 for BAM. Let's go with BAM


#############################################
# The whole kit and kaboodle on a soil sample
import sys
sys.path.append('/nfs1/Koslicki_Lab/koslickd/Repositories/MinHashMetagenomics/src/')
import os
import MinHash as MH
import shutil
import timeit
import h5py
outT0 = timeit.default_timer()
n = 50000
# Get the mins in the training database
fid = h5py.File('/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N'+str(n)+'k31_mins.h5', 'r')
hash_list = set(fid["hash_list"][...])
fid.close()

#Make the CEs for the soil sample
t0 = timeit.default_timer()
soil_sample_file = "/nfs1/Koslicki_Lab/koslickd/CommonKmers/SoilSamples/Data/Metagenomes/4539591.3.fastq"
soil_out_dir = "/nfs1/Koslicki_Lab/koslickd/MinHash/Out/SoilSamples"
MCE_in_comparison = MH.CountEstimator(n=n, max_prime=9999999999971., ksize=31, input_file_name=soil_sample_file, save_kmers='y', hash_list=hash_list)
MCE_in_comparison.export(os.path.join(soil_out_dir, os.path.basename(soil_sample_file)+'_N'+str(n)+'_k31_inComparison.h5'))
MCE_in_all = MH.CountEstimator(n=n, max_prime=9999999999971., ksize=31, input_file_name=soil_sample_file, save_kmers='y')
MCE_in_all.export(os.path.join(soil_out_dir, os.path.basename(soil_sample_file)+'_N'+str(n)+'_k31_all.h5'))
t1 = timeit.default_timer()
print("Sample formation and export time: %f" % (t1-t0))  # 15279.804780 -> 4.24 hours

# Read it back in after it's computed
MCE_in_comparison = MH.import_single_hdf5(os.path.join(soil_out_dir,os.path.basename(soil_sample_file)+'_N'+str(n)+'_k31_inComparison.h5'))
MCE_in_all = MH.import_single_hdf5(os.path.join(soil_out_dir,os.path.basename(soil_sample_file)+'_N'+str(n)+'_k31_all.h5'))

#Read in the training CE's
fid = open('/nfs1/Koslicki_Lab/koslickd/MinHash/Data/FileNames.txt', 'r')
file_names = fid.readlines()
fid.close()
file_names = [name.strip() for name in file_names]
training_n = 5000
out_file_names = ["/nfs1/Koslicki_Lab/koslickd/MinHash/Out/N"+str(training_n)+"k31/" + os.path.basename(item) + ".CE.h5" for item in file_names]
CEs = MH.import_multiple_hdf5(out_file_names)
Y_count_in_comparison = MCE_in_comparison.count_vector(CEs)
eps = .001
(reconstruction, A_eps, A_indicies) = MH.jaccard_count_lsqnonneg(CEs, Y_count_in_comparison, eps)


#reconstruction = reconstruction/float(sum(reconstruction))
#indicies = reconstruction.argsort()[::-1]
#reference_files = []
#for index in indicies:
#    reference_files.append(CEs[index].input_file_name)
#    if reconstruction[index] < eps:  # Add at leat one more index below the threshold (to make sure I don't have an empty list of indexes)
#        break

#print("Number of indicies to build: %d" % len(reference_files))


# Read in the taxonomy and see which are the largest entries of which Y vectors
fid = open('/nfs1/Koslicki_Lab/koslickd/MinHash/Data/Taxonomy.txt', 'r')
taxonomy = fid.readlines()
fid.close()
taxonomy = [item.strip() for item in taxonomy]
taxonomy_names = [item.split('\t')[0] for item in taxonomy]

out_dir = '/scratch/temp/SNAP/training/'
(clusters, LCAs) = MH.cluster_matrix(A_eps, A_indicies, taxonomy, cluster_eps=.01)
training_file_names = MH.make_cluster_fastas(out_dir, LCAs, clusters, CEs)
index_dirs = MH.build_references(training_file_names, out_dir, large_index=True)

#index_dir = "/scratch/temp/SNAP/"
#index_dirs = MH.build_references(reference_files, index_dir)
out_sam = os.path.join(out_dir, "out.sam")
t0 = timeit.default_timer()
MH.stream_aligned_save_unaligned(index_dirs, soil_sample_file, out_sam, format="sam")
t1 = timeit.default_timer()
print("Alignment time: %f" % (t1-t0))

pre, ext = os.path.splitext(out_sam)
out_fastq = pre + ".fastq"
MH.sam2fastq(out_sam, out_fastq)
outT1 = timeit.default_timer()
print("Total time: %f" % (outT1-outT0))


# Let's do it one at a time to see what's going on...(looks like ~100K reads per reference are being filtered out)
MH.top_down_align(soil_sample_file, reference_files, index_dir, "/scratch/temp/SNAP/", save_aligned=False, format="sam")
MH.top_down_align(soil_sample_file, reference_files, index_dir, "/scratch/temp/SNAP/", save_aligned=True, format="bam")








# Clean up
os.remove(out_sam)
for index_dir in index_dirs:
    shutil.rmtree(index_dir)













#############
# SAM/BAM ugliness
#/local/cluster/samtools/bin/./samtools view -b -f 4 out.sam | more

#/local/cluster/jdk1.8.0_71/bin/java -jar picard.jar

#Will need to look into the details for paired reads, as I want to count as aligned any pair of reads where one of the pairs mapped (OR DO I?!)
# This is giving me error...ugh! Just run snap twice with different -F values
#/local/cluster/jdk1.8.0_71/bin/java -jar /local/cluster/picard/dist/picard.jar FilterSamReads I=/nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/SNAP/out.bam O=/nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/SNAP/out_aligned.bam FILTER=includeAligned
#/local/cluster/jdk1.8.0_71/bin/java -jar /local/cluster/picard/dist/picard.jar FilterSamReads I=/nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/SNAP/out.bam O=/nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/SNAP/out_unaligned.bam FILTER=excludeAligned

# all
# samtools view 4539591.3.fastq_G000502875.fna_aligned.bam | wc -l
# 6033053

# read unmapped
# samtools view -f4 4539591.3.fastq_G000502875.fna_aligned.bam | wc -l
# 6000248

# read mapped
# samtools view -F4 4539591.3.fastq_G000502875.fna_aligned.bam | wc -l
# 32805

# view alignment
# samtools sort -@ 20 4539591.3.fastq_G000502875.fna_aligned.bam -o test.bam
# samtools index test.bam
# samtools faidx G000502875.fna
# samtools tview test.bam

# Let's try to see why it's showing so many reads unmapped
# snap-aligner single G000502875 /nfs1/Koslicki_Lab/koslickd/CommonKmers/SoilSamples/Data/Metagenomes/4539591.3.fastq -o test.bam -F a
# Total Reads    Aligned, MAPQ >= 10    Aligned, MAPQ < 10     Unaligned              Too Short/Too Many Ns  Filtered                  Reads/s   Time in Aligner (s)
# 65,927,593     114,525 (0.17%)        2,029 (0.00%)          0 (0.00%)              6,000,248 (9.10%)      59,810,791 (90.72%)       688,668   96
# samtools view test.bam | wc -l
# 6116802 == 114525+2029+6000248
# samtools view -f4 test.bam | wc -l
# 6000248
# samtools view -F4 test.bam | wc -l
# 116554 == 114525+2029





#snap-aligner single G000502875 /nfs1/Koslicki_Lab/koslickd/CommonKmers/SoilSamples/Data/Metagenomes/4539591.3.fastq -o test.bam -F a -f
# Total Reads    Aligned, MAPQ >= 10    Aligned, MAPQ < 10     Unaligned              Too Short/Too Many Ns  Filtered                  Reads/s   Time in Aligner (s)
# 65,927,593     0 (0.00%)              116,554 (0.18%)        0 (0.00%)              6,000,248 (9.10%)      59,810,791 (90.72%)       691,884   95

#snap-aligner single G000502875 /nfs1/Koslicki_Lab/koslickd/CommonKmers/SoilSamples/Data/Metagenomes/4539591.3.fastq -o test.bam -F a -f -d 14 -mrl 50
# Total Reads    Aligned, MAPQ >= 10    Aligned, MAPQ < 10     Unaligned              Too Short/Too Many Ns  Filtered                  Reads/s   Time in Aligner (s)
# 65,927,593     0 (0.00%)              116,554 (0.18%)        0 (0.00%)              6,000,248 (9.10%)      59,810,791 (90.72%)       691,536   95

# snap-aligner single G000502875 /nfs1/Koslicki_Lab/koslickd/CommonKmers/SoilSamples/Data/Metagenomes/4539591.3.fastq -o test.bam -E sm -f
# samtools view -f4 test.bam | wc -l
# 116554


# sort reads
# cat 4539591.3.fastq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > sorted.fastq

# Paired: save the pairs where at least one of the pairs aligned
# samtools merge merged.bam <(samtools view -hb -F4 -f8 test.sam) <(samtools view -hb -F8 -f4 test.sam) <(samtools view -hb -F12 test.sam) -f
# Paired: save the pairs where neither one of the pairs aligned
# samtools merge merged.bam <(samtools view -hb -f12 test.sam) -f

# Better yet!
# samtools view -hb -f12 test.sam -o unaligned.sam -U aligned.sam




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




###################
# Taxonomy
(clusters, LCAs) = MH.cluster_matrix(A_eps, A_indicies, taxonomy, cluster_eps=.01)


###################
# Form the training data
import screed
out_dir = '/scratch/temp/SNAP/training/'
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

out_file_names = []
for cluster_index in range(len(clusters)):
    cluster = clusters[cluster_index]
    LCA = LCAs[cluster_index]
    out_file_name = os.path.join(out_dir, LCA + "_" + str(cluster_index) + "_" + ".fa")  # put the cluster index in the name in case there are shared LCAs
    out_file_names.append(out_file_name)
    out_file = open(out_file_name, 'w')
    for index in cluster:
        file_name = CEs[index].input_file_name
        for record in screed.open(file_name):
            out_file.write(">" + LCA)
            out_file.write("\n")
            out_file.write(record.sequence)
            out_file.write("\n")
    out_file.close()

def _write_single(tup):
    out_dir = tup[0]
    LCA = tup[1]
    cluster_index = tup[2]
    input_file_names = tup[3]
    out_file_name = os.path.join(out_dir, LCA + "_" + str(cluster_index) + "_" + ".fa")  # put the cluster index in the name in case there are shared LCAs
    out_file = open(out_file_name, 'w')
    for file_name in input_file_names:
        for record in screed.open(file_name):
            out_file.write(">" + LCA)
            out_file.write("\n")
            out_file.write(record.sequence)
            out_file.write("\n")
    out_file.close()
    return out_file_name

write_single((out_dir, LCAs[0], 0, [CEs[i].input_file_name for i in clusters[0]]))


pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
res = pool.map(_write_single, zip(repeat(out_dir), LCAs, range(len(LCAs)), [[CEs[i].input_file_name for i in cluster] for cluster in clusters]), chunksize=1)



for record in screed.open(self.input_file_name):
            self.add_sequence(record.sequence)













# Export the clustered matrix
clustered_inds = []
for cluster in clusters:
    for ind in cluster:
        clustered_inds.append(ind)

A_clustered = A_eps[clustered_inds,:][:,clustered_inds]


############
# All index build time: 70.36 min. (1227 bacterial genomes)
# All align time: 235s (snap-aligner single . ../4539591.3.fastq -o test_all.bam)
# Total Reads    Aligned, MAPQ >= 10    Aligned, MAPQ < 10     Unaligned              Too Short/Too Many Ns     Reads/s   Time in Aligner (s)
# 65,927,593     3,078,885 (4.67%)      1,562,255 (2.37%)      55,286,205 (83.86%)    6,000,248 (9.10%)         279,994   235
# All align time with multiple hits:  (snap-aligner single . ../4539591.3.fastq -o test_all.bam -om 5 -D 5 -d 20)
# Total Reads    Aligned, MAPQ >= 10    Aligned, MAPQ < 10     Unaligned              Too Short/Too Many Ns   Extra Alignments    Reads/s   Time in Aligner (s)
# 65,927,593     6,544,469 (9.93%)      4,004,920 (6.07%)      49,377,956 (74.90%)    6,000,248 (9.10%)       99,489,436          209,932   314