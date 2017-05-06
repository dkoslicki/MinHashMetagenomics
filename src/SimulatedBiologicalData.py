# This script will simulate some metagenomic data and compare:
# 1. Brute force ground truth
# 2. Classic min hash
# 3. Containment min hash
# via Gem-Sim
import os
import screed
import numpy as np
import subprocess
from pybloom import BloomFilter  # basic bloom filter for comparison purposes
import khmer
import MinHash as MH
import bz2
import tempfile
import matplotlib.pyplot as plt
from multiprocessing.dummy import Pool

num_threads = 8
num_genomes = 20
num_reads = 1000000
num_replicates = 20
python_loc = "python"
gen_sim_loc = "/home/dkoslicki/Documents/GemSIM_v1.6/GemReads.py"
prime = 9999999999971  # taking hashes mod this prime
p = 0.001  # false positive rate for the bloom filter
ksize = 11  # k-mer length
max_h = 5000  # Maximum number of hashes to use
hash_range = range(max_h, 10, -100)  # range of hash values to test


def make_simulation(num_genomes, num_reads, python_loc, gen_sim_loc):
	# First, read in the available files
	file_names = list()
	fid = open(os.path.abspath('../data/Genomes/FileNames.txt'), 'r')
	for line in fid.readlines():
		file_names.append(os.path.abspath(os.path.join('../data/Genomes/', os.path.basename(line.strip()))))

	fid.close()

	# Make a random selection of the genomes
	selected_genomes = np.random.choice(file_names, num_genomes, replace=False)

	# Make an abundance file with random abundances, use uniform distribution for now
	abundances = np.random.rand(num_genomes)
	abundances /= np.sum(abundances)
	abundances_temp_file = tempfile.NamedTemporaryFile('w', suffix=".txt")
	abundances_temp_file_name = abundances_temp_file.name
	abundances_temp_file.close()
	fid = open(abundances_temp_file_name, 'w')
	for i in range(num_genomes):
		fid.write('%s\t%f\n' % (os.path.basename(selected_genomes[i]), abundances[i]))

	fid.close()

	# Run Gen Sim
	genomes_loc = '../data/Genomes/'
	cmd = python_loc + " " + gen_sim_loc + " -R " + genomes_loc + " -a " + abundances_temp_file_name + " -n " + str(num_reads) + \
		" -l 100 -c -q 64 -o " + os.path.join(os.path.dirname(abundances_temp_file_name), os.path.splitext(os.path.basename(abundances_temp_file_name))[0]) + " -m " + \
		os.path.join(os.path.dirname(gen_sim_loc), "models/ill100v5_s.gzip")

	devnull = open(os.devnull, 'w')  # Gen sim has a "divide by zero" warning, so let's suppress it
	test = subprocess.check_output(cmd, shell=True, stderr=devnull)
	return os.path.join(os.path.dirname(abundances_temp_file_name), os.path.splitext(os.path.basename(abundances_temp_file_name))[0]) + "_single.fastq", abundances_temp_file_name, selected_genomes


def create_relative_errors(num_genomes, num_reads, python_loc, gen_sim_loc, prime, p, ksize, hash_range):
	# Make a simulation
	simulation_file, abundances_file, selected_genomes = make_simulation(num_genomes, num_reads, python_loc, gen_sim_loc)

	# Get simulation k-mers, use canonical k-mers
	# Simultaneously, make the min hash sketch of the simulation
	simulation_kmers = set()
	simulation_MHS = MH.CountEstimator(n=max_h, max_prime=prime, ksize=ksize, save_kmers='y')
	for record in screed.open(simulation_file):
		seq = record.sequence
		for i in xrange(len(seq) - ksize + 1):
			kmer = seq[i:i+ksize]
			kmer_rev = khmer.reverse_complement(kmer)
			if kmer < kmer_rev:
				simulation_kmers.add(kmer)
				simulation_MHS.add(kmer)
			else:
				simulation_kmers.add(kmer_rev)
				simulation_MHS.add(kmer_rev)

	# Use them to populate a bloom filter
	simulation_bloom = BloomFilter(capacity=1.1*len(simulation_kmers), error_rate=p)
	simulation_kmers_length = len(simulation_kmers)  # in practice, this would be computed when the bloom filter is created
	# or can use an estimate based on the bloom filter entries
	for kmer in simulation_kmers:
		simulation_bloom.add(kmer)

	# Use pre-computed data to load the kmers and the sketches
	base_names = [os.path.basename(item) for item in selected_genomes]
	# Load the sketches
	genome_sketches = MH.import_multiple_from_single_hdf5(os.path.abspath('../data/Genomes/AllSketches.h5'), base_names)
	# Get the true number of kmers
	genome_lengths = list()
	for i in range(len(genome_sketches)):
		genome_lengths.append(genome_sketches[i]._true_num_kmers)

	# Get *all* the kmers for computation of ground truth
	genome_kmers = list()
	for i in range(len(base_names)):
		name = base_names[i]
		kmers = set()
		fid = bz2.BZ2File(os.path.abspath(os.path.join('../data/Genomes/', name + ".kmers.bz2")), 'r')
		for line in fid.readlines():
			kmers.add(line.strip())
		fid.close()
		genome_kmers.append(kmers)

	# Calculate the true Jaccard index
	true_jaccards = list()
	for kmers in genome_kmers:
		true_jaccard = len(kmers.intersection(simulation_kmers)) / float(len(kmers.union(simulation_kmers)))
		true_jaccards.append(true_jaccard)

	# Calculate the min hash estimate of jaccard index
	MH_relative_errors = list()
	CMH_relative_errors = list()
	for h in hash_range:
		MH_jaccards = list()
		for MHS in genome_sketches:
			# Down sample each sketch to h
			MHS.down_sample(h)
			simulation_MHS.down_sample(h)
			MH_jaccard = MHS.jaccard(simulation_MHS)
			MH_jaccards.append(MH_jaccard)

		MH_jaccards_corrected = list()
		for MHS in genome_sketches:
			MHS_set = set(MHS._mins)
			sample_set = set(simulation_MHS._mins)
			MH_jaccard = len(set(list(MHS_set.union(sample_set))[0:h]).intersection(MHS_set.intersection(sample_set))) / float(h)
			MH_jaccards_corrected.append(MH_jaccard)

		# Calculate the containment min hash estimate of the jaccard index
		CMH_jaccards = list()
		for i in range(len(genome_sketches)):
			genome_kmers_len = genome_lengths[i]  # pre-computed when creating the "training" data
			MHS = genome_sketches[i]
			# down sample each sketch to h
			MHS.down_sample(h)
			kmers = MHS._kmers  # use only the k-mers in the min hash sketch
			int_est = 0
			for kmer in kmers:
				if kmer in simulation_bloom:  # test if the k-mers are in the simulation bloom filter
					int_est += 1
			int_est -= p*h  # adjust for false positive rate
			containment_est = int_est / float(h)
			containment_est_jaccard = genome_kmers_len * containment_est / \
				(genome_kmers_len + simulation_kmers_length - genome_kmers_len * containment_est)
			CMH_jaccards.append(containment_est_jaccard)

		# compute the average deviation from the truth (relative error)
		true_jaccards = np.array(true_jaccards)
		MH_jaccards = np.array(MH_jaccards)
		CMH_jaccards = np.array(CMH_jaccards)
		MH_mean = np.mean(np.abs(true_jaccards - MH_jaccards)/true_jaccards)
		CMH_mean = np.mean(np.abs(true_jaccards - CMH_jaccards)/true_jaccards)
		#print("Classic min hash mean relative error: %f" % MH_mean)
		#print("Containment min hash mean relative error: %f" % CMH_mean)
		MH_relative_errors.append(MH_mean)
		CMH_relative_errors.append(CMH_mean)

	# remove temp files
	os.remove(simulation_file)
	os.remove(abundances_file)
	# return the relative errors
	return MH_relative_errors, CMH_relative_errors, simulation_kmers_length, np.mean(genome_lengths)


# Now let's parallelize it and do a bunch of replicates
def dummy_wrapper(arg):
	return create_relative_errors(*arg)

pool = Pool(processes=num_threads)
res = pool.map(dummy_wrapper, [(num_genomes, num_reads, python_loc, gen_sim_loc, prime, p, ksize, hash_range)
								for item in range(num_replicates)])

# Plot the results
MH_results = np.zeros((num_replicates, len(hash_range)))
CMH_results = np.zeros((num_replicates, len(hash_range)))
simulation_kmers_lengths = np.zeros(num_replicates)
mean_genome_lengths = np.zeros(num_replicates)
for i in range(num_replicates):
	MH_results[i, :] = res[i][0]
	CMH_results[i, :] = res[i][1]
	simulation_kmers_lengths[i] = res[i][2]
	mean_genome_lengths[i] = res[i][3]

np.savetxt(os.path.abspath('../data/SimulatedMetagenomes/MinHash_results.txt'), MH_results)
np.savetxt(os.path.abspath('../data/SimulatedMetagenomes/ContainmentMinHash_results.txt'), CMH_results)

font = {'family': 'serif',
		'color':  'black',
		'weight': 'normal',
		'size': 18,
		}

plt.figure()
plt.errorbar(hash_range, np.mean(MH_results, 0), yerr=np.std(MH_results, 0), fmt='--r', ecolor=[1, 0, 0, .2], label="Classic Min Hash")
plt.errorbar(hash_range, np.mean(CMH_results, 0), yerr=np.std(CMH_results, 0), fmt='b', ecolor=[0, 0, 1, .2], label="Containment Min Hash")
axes = plt.gca()
axes.text(-.1, 1, 'b)', horizontalalignment='left', verticalalignment='bottom', fontdict=font, transform=axes.transAxes)
plt.ylabel('Relative error')
plt.xlabel('Number of hashes')
plt.xlim([min(hash_range), max(hash_range)])
ticks = plt.xticks()
ticks[0][0] = min(hash_range)
plt.xticks(ticks[0])
plt.legend()

plt.savefig('../Paper/Figs/SimulatedBiologicalData.png')

# Save parameters for input to LaTeX paper
fid = open(os.path.abspath('../Paper/Data/SimulatedBiologicalDataNumGenomes.txt'), 'w')
fid.write('%d' % num_genomes)
fid.close()
fid = open(os.path.abspath('../Paper/Data/SimulatedBiologicalDataNumReads.txt'), 'w')
fid.write('%d' % num_reads)
fid.close()
fid = open(os.path.abspath('../Paper/Data/SimulatedBiologicalDataNumReplicates.txt'), 'w')
fid.write('%d' % num_replicates)
fid.close()
fid = open(os.path.abspath('../Paper/Data/SimulatedBiologicalDatap.txt'), 'w')
fid.write('%f' % p)
fid.close()
fid = open(os.path.abspath('../Paper/Data/SimulatedBiologicalDataksize.txt'), 'w')
fid.write('%d' % ksize)
fid.close()
fid = open(os.path.abspath('../Paper/Data/SimulatedBiologicalData_rel_size.txt'), 'w')
fid.write('%3.3f\%%' % float(100*np.mean(simulation_kmers_lengths) / float(np.mean(mean_genome_lengths))))
fid.close()




