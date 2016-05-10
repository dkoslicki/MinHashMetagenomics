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


def build_reference(reference_file, output_dir):
    cmd = "/home/pi/koslickd/./snap-aligner index " + reference_file + " " + output_dir
    dummy = subprocess.call(cmd, shell=True,  stdout=FNULL, stderr=subprocess.STDOUT)
    return


def align_reads(index_dir, sample_file, out_file, aligned=True):  # NOTE: snap-aligner will take SAM and BAM as INPUT!!
    if aligned:
        cmd = "/home/pi/koslickd/./snap-aligner single " + index_dir + " " + sample_file + " -o " + out_file + " -F a -t 48 -d 20"
    elif not aligned:
        cmd = "/home/pi/koslickd/./snap-aligner single " + index_dir + " " + sample_file + " -o " + out_file + " -F u -t 48 -d 20"
    else:
        raise Exception("aligned must be True or False")
    dummy = subprocess.call(cmd, shell=True,  stdout=FNULL, stderr=subprocess.STDOUT)
    return


def sam2fastq(sam_file, out_file):
    cmd = "cat " + sam_file + " | grep -v ^@ | awk '{print \"@\"$1\"\\n\"$10\"\\n+\\n\"$11}' > " + out_file
    dummy = subprocess.call(cmd, shell=True,  stdout=FNULL, stderr=subprocess.STDOUT)
    return

build_reference(reference_file, index_dir)
align_reads(index_dir, sample_file, sam_out_file, aligned=True)
sam2fastq(sam_out_file, fastq_out_basename + "_aligned.fastq")
align_reads(index_dir, sample_file, sam_out_file, aligned=False)
sam2fastq(sam_out_file, fastq_out_basename + "_unaligned.fastq")
os.remove(sam_out_file)

# More efficient thing to do is figure out how to get the output SAM to split over aligned and unaligned, rather than reading it twice...
# Since snap-aligner also accepts streaming inputs, if I didn't care to save the aligned reads, I could just chain a bunch of calls to it together with different references!

# Do a bunch of them
reference_files = []
reconstruction = MH.jaccard_count_lsqnonneg(CEs, Y_count_in_comparison, eps)
N = 10
indicies = reconstruction.argsort()[-N:][::-1]
for index in indicies:
    reference_files.append(CEs[index].input_file_name)

reference_files = ['/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000008565.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000015425.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000012825.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000006885.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000016965.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000008525.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000012905.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000387985.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000008805.fna', '/nfs1/Koslicki_Lab/koslickd/CommonKmers/TrainingOnRepoPhlAnPython/TrainDataIn/Genomes/G000154225.fna']
for i in range(len(reference_files)):
    reference_file = reference_files[i]
    if i == 0:
        fastq_out_basename = os.path.join(out_dir, os.path.basename(sample_file) + "_" + os.path.basename(reference_file))
        build_reference(reference_file, index_dir)
        align_reads(index_dir, sample_file, sam_out_file, aligned=True)
        sam2fastq(sam_out_file, fastq_out_basename + "_aligned.fastq")
        align_reads(index_dir, sample_file, sam_out_file, aligned=False)
        sam2fastq(sam_out_file, fastq_out_basename + "_unaligned.fastq")  # Better yet, just call this unaligned.fastq
        os.remove(sam_out_file)
    else:
        fastq_out_basename = os.path.join(out_dir, os.path.basename(sample_file) + "_" + os.path.basename(reference_file))
        prev_fastq_file = os.path.join(out_dir, os.path.basename(sample_file) + "_" + os.path.basename(reference_files[i-1]) + "_unaligned.fastq")  # Better yet, just call this unaligned.fastq
        build_reference(reference_file, index_dir)
        align_reads(index_dir, prev_fastq_file, sam_out_file, aligned=True)
        sam2fastq(sam_out_file, fastq_out_basename + "_aligned.fastq")
        align_reads(index_dir, prev_fastq_file, sam_out_file, aligned=False)
        sam2fastq(sam_out_file, fastq_out_basename + "_unaligned.fastq")
        os.remove(sam_out_file)
        os.remove(prev_fastq_file)




#############
# SAM/BAM ugliness
#/local/cluster/samtools/bin/./samtools view -b -f 4 out.sam | more

#/local/cluster/jdk1.8.0_71/bin/java -jar picard.jar

#Will need to look into the details for paired reads, as I want to count as aligned any pair of reads where one of the pairs mapped (OR DO I?!)
# This is giving me error...ugh! Just run snap twice with different -F values
#/local/cluster/jdk1.8.0_71/bin/java -jar /local/cluster/picard/dist/picard.jar FilterSamReads I=/nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/SNAP/out.bam O=/nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/SNAP/out_aligned.bam FILTER=includeAligned
#/local/cluster/jdk1.8.0_71/bin/java -jar /local/cluster/picard/dist/picard.jar FilterSamReads I=/nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/SNAP/out.bam O=/nfs1/Koslicki_Lab/koslickd/MinHash/Out/Temp/SNAP/out_unaligned.bam FILTER=excludeAligned
















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
