# This script will query the sketched k-mers against a soil metagenome bloom filter
import subprocess
import os
import numpy as np
import MinHash as MH

ksize = 21  # k-mer length
max_h = 500  # max number of hashes in sketch
query_per_sequence_loc = 'QueryPerSequence/./query_per_sequence'
metagenome_bloom_filter = '/home/dkoslicki/Dropbox/Repositories/MinHash/data/4539585.3.fastq.21mers.bc.p.01'
metagenome_kmers_total = 916485607  # via jellyfish stats on a count bloom filter (need to streamline this)


# Get virus names
file_names = list()
fid = open(os.path.abspath('../data/Viruses/FileNames.txt'), 'r')
for line in fid.readlines():
	file_names.append(os.path.abspath(os.path.join('../data/Viruses/', os.path.basename(line.strip()))))

fid.close()

# Get all the hashes so we know the size
base_names = [os.path.basename(item) for item in file_names]
genome_sketches = MH.import_multiple_from_single_hdf5(os.path.abspath('../data/Viruses/AllSketches.h5'), base_names)

# query the bloom filter
jaccards = np.zeros(len(file_names))
it = 0
for genome in file_names:
	name = os.path.basename(genome)
	CMH = genome_sketches[it]
	genome_kmers_len = CMH._true_num_kmers
	cmd = query_per_sequence_loc + " " + metagenome_bloom_filter + " " + \
		os.path.abspath(os.path.join('../data/Viruses/', name + ".Hash" + str(ksize) + "mers.fa"))
	int_est = int(subprocess.check_output(cmd, shell=True))
	int_est -= 0.01*int_est
	containment_est = int_est / float(max_h)
	containment_est_jaccard = genome_kmers_len * containment_est / \
							(genome_kmers_len + metagenome_kmers_total - genome_kmers_len * containment_est)
	jaccards[it] = containment_est_jaccard
	if it % 100 == 0:
		print('On iteration: %d' % it)
	it += 1

np.savetxt(os.path.abspath('../data/4539585.3.fastq.jaccards.txt'), jaccards)


## What was the largest containment guy
jaccards = np.loadtxt(os.path.abspath('../data/4539585.3.fastq.jaccards.txt'))
it = 0
containments = list()
for genome in file_names:
	name = os.path.basename(genome)
	CMH = genome_sketches[it]
	genome_kmers_len = CMH._true_num_kmers
	j = jaccards[it]
	g = genome_kmers_len
	m = metagenome_kmers_total
	containment_est = (j*(g+m))/float(g*(1+j))
	it += 1
	containments.append(containment_est*max_h)

