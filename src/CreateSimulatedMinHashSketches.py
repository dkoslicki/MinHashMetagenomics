# Here we will pre-compute the min hash sketches to speed the simulation. Could parallelize this
# to make it go faster, but it only needs to be done once, so I don't really care
import os
import screed
import khmer
import MinHash as MH
import bz2

prime = 9999999999971  # taking hashes mod this prime
ksize = 11  # k-mer length
max_h = 5000  # max number of hashes in sketch

file_names = list()
fid = open(os.path.abspath('../data/Genomes/FileNames.txt'), 'r')
for line in fid.readlines():
	file_names.append(os.path.abspath(os.path.join('../data/Genomes/', os.path.basename(line.strip()))))

fid.close()

# Get genome k-mers and also make min hash sketches
genome_sketches = list()
for genome in file_names:
	kmers = set()
	name = os.path.basename(genome)
	MHS = MH.CountEstimator(n=max_h, max_prime=prime, ksize=ksize, save_kmers='y')
	for record in screed.open(genome):
		seq = record.sequence
		for i in xrange(len(seq) - ksize + 1):
			kmer = seq[i:i+ksize]
			kmer_rev = khmer.reverse_complement(kmer)
			if kmer < kmer_rev:
				kmers.add(kmer)
				MHS.add(kmer)
			else:
				kmers.add(kmer_rev)
				MHS.add(kmer_rev)
	MHS._true_num_kmers = len(kmers)
	MHS.input_file_name = os.path.basename(genome)
	genome_sketches.append(MHS)
	# export the kmers
	fid = bz2.BZ2File(os.path.abspath(os.path.join('../data/Genomes/', name + ".kmers.bz2")), 'w')
	for kmer in kmers:
		fid.write("%s\n" % kmer)
	fid.close()

# Export all the sketches
base_names = [os.path.basename(item) for item in file_names]
MH.export_multiple_to_single_hdf5(genome_sketches, os.path.abspath('../data/Genomes/AllSketches.h5'))
