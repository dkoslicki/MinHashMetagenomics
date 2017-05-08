# This script will query the sketched k-mers against a soil metagenome bloom filter
import subprocess
import os
import numpy as np
import MinHash as MH
from itertools import *
from multiprocessing.dummy import Pool
import multiprocessing
import matplotlib.pyplot as plt

num_threads = multiprocessing.cpu_count()
ksize = 21  # k-mer length
max_h = 500  # max number of hashes in sketch
p = 0.01
query_per_sequence_loc = os.path.abspath('QueryPerSequence/./query_per_sequence')
metagenome_bloom_filter = os.path.abspath('../data/MetagenomeBloom.jf')
#metagenome_kmers_total = 916485607  # via jellyfish stats on a count bloom filter (need to streamline this)
# Note that this number changed since I am restricting myself to the paired reads (no orphaned guys)
fid = open(os.path.abspath("../Paper/Data/MetagenomeTotalKmers.txt"), 'r')
metagenome_kmers_total = int(fid.readlines()[0].strip())
fid.close()

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
def count_jaccard(genome_sketch, file_name, ksize, p, max_h, metagenome_kmers_total):
	name = os.path.basename(file_name)
	CMH = genome_sketch
	genome_kmers_len = CMH._true_num_kmers
	cmd = query_per_sequence_loc + " " + metagenome_bloom_filter + " " + \
		os.path.abspath(os.path.join('../data/Viruses/', name + ".Hash" + str(ksize) + "mers.fa"))
	int_est = int(subprocess.check_output(cmd, shell=True))
	int_est -= p*int_est
	containment_est = int_est / float(max_h)
	containment_est_jaccard = genome_kmers_len * containment_est / \
							(genome_kmers_len + metagenome_kmers_total - genome_kmers_len * containment_est)
	return containment_est_jaccard

def unwrap_make_jaccard(arg):
	return count_jaccard(*arg)

pool = Pool(processes=num_threads)
jaccards = np.array(pool.map(unwrap_make_jaccard, zip(genome_sketches, file_names, repeat(ksize), repeat(p), repeat(max_h), repeat(metagenome_kmers_total))))

np.savetxt(os.path.abspath('../data/MetagenomeViruses.jaccards.txt'), jaccards)


# What was the largest containment guy (If I want an estimate of the pure cardinality of the intersection,
# then use containments.append(containment_est*max_h)
it = 0
containments = list()
for genome in file_names:
	CMH = genome_sketches[it]
	genome_kmers_len = CMH._true_num_kmers
	j = jaccards[it]
	g = genome_kmers_len
	m = metagenome_kmers_total
	containment_est = (j*(g+m))/float(g*(1+j))
	it += 1
	containments.append(containment_est)

#plt.figure()
#plt.plot(containments)
pos = np.array(containments).argmax()
name = os.path.basename(file_names[pos])
fid = open(file_names[pos], 'r')
head = fid.readline()
fid.close()
accession = head.split(' ')[0][1:]
head = ' '.join(head.split(',')[0].split(' ')[1:])
fid = open(os.path.abspath('../Paper/Data/FoundOrganismName.txt'), 'w')
fid.write("%s" % head)
fid.close()
fid = open(os.path.abspath('../Paper/Data/FoundOrganismAccession.txt'), 'w')
fid.write("%s" % accession)
fid.close()
fid = open(os.path.abspath('../Paper/Data/FoundOrganismFileName.txt'), 'w')
fid.write("%s" % file_names[pos])
fid.close()
fid = open(os.path.abspath('../Paper/Data/FoundOrganismContainment.txt'), 'w')
fid.write("%1.4f" % containments[pos])
fid.close()
fid = open(os.path.abspath('../Paper/Data/FoundOrganismJaccard.txt'), 'w')
fid.write("%.3e" % jaccards[pos])
fid.close()


# What was the largest jaccard, don't think this is a really good way to get the guys in there
# coverage is a better estimate of presence/absence (jaccard just returns megaviruses)
#pos = jaccards.argmax()
#name = os.path.basename(file_names[pos])
#plt.figure()
#plt.plot(jaccards)
