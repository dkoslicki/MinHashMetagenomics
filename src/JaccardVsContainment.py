# This will be a script to generate a figure demonstrating the advantage of the containment approach vs
# the typical min hash jaccard approach when the set sizes are very different.
import MinHash as MH  # Our implementation of Min Hash
import numpy as np
import matplotlib.pyplot as plt
#from pybloom import BloomFilter  # basic bloom filter for comparison purposes (only for Python2)
from pybloom_live import BloomFilter  # basic bloom filter for comparison purposes
import os

prime = 9999999999971  # taking hashes mod this prime
p = 0.001  # false positive rate for the bloom filter
n1 = 10000  # seq1 sequence length
n2 = 15  # seq2 sequence length
ksize = 11  # k-mer length
h = 100  # number of hashes in sketch


#  function to return even range of jaccard values, j is desired jaccard, returns length of string you need to append
def i_size(j, k, n):
	return int(n*j/float(1-j) + k)

i_range = [i_size(val, ksize, n1) for val in np.arange(0, 1, 0.005)]
#i_range = [i_size(val, ksize, n1) for val in np.arange(0, 1, 0.05)]  # For quick debugging
true_jaccards = np.zeros(len(i_range))
estimate_jaccards = np.zeros(len(i_range))
containment_jaccards = np.zeros(len(i_range))
it = 0

for i_size in i_range:
	# Append a common string to two different random strings (true jaccard will be ~ i_size/n)
	common_string = ''.join(np.random.choice(['A', 'C', 'T', 'G'], i_size))
	seq1 = ''.join(np.random.choice(['A', 'C', 'T', 'G'], n1)) + common_string
	# Make seq2 a smaller sequence than seq1
	seq2 = ''.join(np.random.choice(['A', 'C', 'T', 'G'], n2)) + common_string

	# Calculate exact Jaccard index
	kmers1 = set()
	kmers2 = set()
	for i in range(len(seq1) - ksize + 1):
		kmers1.add(seq1[i:i+ksize])

	for i in range(len(seq2) - ksize + 1):
		kmers2.add(seq2[i:i+ksize])

	true_jaccard = len(kmers1.intersection(kmers2)) / float(len(kmers1.union(kmers2)))
	true_jaccards[it] = true_jaccard

	# Calculate sourmash estimate of Jaccard index
	E1 = MH.CountEstimator(n=h, max_prime=prime, ksize=ksize, save_kmers='y')
	E2 = MH.CountEstimator(n=h, max_prime=prime, ksize=ksize, save_kmers='y')
	E1.add_sequence(seq1)
	E2.add_sequence(seq2)
	estimate_jaccard = E1.jaccard(E2)
	estimate_jaccards[it] = estimate_jaccard

	# Containment version.
	# Bloom filter
	f = BloomFilter(capacity=i_size+n1, error_rate=p)
	len_kmers_1 = 0
	for val in kmers1:
		if val not in f:
			len_kmers_1 += 1
			f.add(val)
	#len_kmers_1 *= (1 - p)  # adjust for the false positive rate, shouldn't need to do this as I'm just adding elements
	int_est = 0
	for val in E2._kmers:
		#if val in f:  # in python2, no distinguishing between byte and utf-8 string
		if val is not '':
			if val.decode("utf-8") in f:
				int_est += 1
	int_est -= p*h  # adjust for the false positive rate
	containment_est = int_est / float(h)

	# Calculate the containment estimate of jaccard, len(kmers2) is exact (as in practice this is part of the training
	# database and so only needs to be done once (and the genomes are relatively small so this is no big cost)
	containment_est_jaccard = \
		len(kmers2) * containment_est / \
		(len(kmers2) + len_kmers_1 - len(kmers2) * containment_est)

	containment_jaccards[it] = containment_est_jaccard
	it += 1

font = {'family': 'serif',
		'color':  'black',
		'weight': 'normal',
		'size': 18,
		}

differences = true_jaccards - estimate_jaccards
sorted_true = sorted(true_jaccards)
sorted_estimates = np.array([x for (y, x) in sorted(zip(true_jaccards, estimate_jaccards), key=lambda pair: pair[0])])
sorted_differences = sorted_true - sorted_estimates
plt.figure()
plt.plot(sorted_true, sorted_differences)
axes = plt.gca()
axes.set_ylim([np.min(plt.yticks()[0])*1.5, np.max(plt.yticks()[0])*1.5])
plt.title('True - estimate Jaccard index')
plt.text(0, 0, 'Underestimate', rotation=90, horizontalalignment='center', verticalalignment='bottom', multialignment='center', color='b', fontsize=14)
plt.text(0, 0, 'Overestimate', rotation=90, horizontalalignment='center', verticalalignment='top', multialignment='center', color='r', fontsize=14)
plt.axhline(0, color='black', linestyle='dashed', linewidth=2)
plt.ylabel('Difference')
plt.xlabel('True Jaccard index')
plt.savefig('../Paper/Figs/Differences.png')

# Do a true vs sourmash estimate plot
plt.figure()
f, ax = plt.subplots()
ax.plot([0, 1], [0, 1], ls="--", c=".3")
ax.plot(sorted_true, sorted_estimates)
plt.ylabel('Estimate Jaccard')
plt.xlabel('True Jaccard')
plt.title('Classic Min Hash')
axes = plt.gca()
#axes.text(-.2, 1.15, 'a)', horizontalalignment='left', verticalalignment='bottom', fontdict=font)
plt.savefig('../Paper/Figs/TrueVsEstimate.png')

# Do a relative error plot
plt.figure()
plt.plot(sorted_true, sorted_differences / sorted_true)
axes = plt.gca()
axes.set_ylim([-1, 1])
plt.axhline(0, color='black', linestyle='dashed', linewidth=2)
plt.ylabel('Relative error')
plt.xlabel('True Jaccard index')
plt.savefig('../Paper/Figs/RelativeError.png')

plt.figure()
n, bins, patches = plt.hist(differences, 50, normed=1, facecolor='green', alpha=0.75)
plt.axvline(0, color='b', linestyle='dashed', linewidth=2)
plt.axvline(np.mean(differences), color='black', linestyle='dashed', linewidth=2)
plt.title('Histogram of (true - estimate) Jaccard index\n Mean: %f' % np.mean(differences))
plt.text(0, max(plt.yticks()[0])-1, 'Underestimate', rotation=0, horizontalalignment='left', verticalalignment='top', multialignment='left', color='b', fontsize=14)
plt.text(plt.xticks()[0][1], max(plt.yticks()[0])-1, 'Overestimate', rotation=0, horizontalalignment='left', verticalalignment='top', multialignment='left', color='r', fontsize=14)
plt.xlabel('Difference')
plt.savefig('../Paper/Figs/Histogram.png')

# Containment guys
plt.figure()
sorted_true = sorted(true_jaccards)
sorted_containment_estimates = np.array([x for (y, x) in sorted(zip(true_jaccards, containment_jaccards), key=lambda pair: pair[0])])
containment_differences = sorted_true - sorted_containment_estimates
n, bins, patches = plt.hist(containment_differences, 50, normed=1, facecolor='green', alpha=0.75)
plt.axvline(0, color='b', linestyle='dashed', linewidth=2)
plt.axvline(np.mean(containment_differences), color='black', linestyle='dashed', linewidth=2)
plt.title('Histogram of (true - corrected estimate) Jaccard index\n Mean: %f' % np.mean(containment_differences))
plt.text(0, max(plt.yticks()[0])-1, 'Underestimate', rotation=0, horizontalalignment='left', verticalalignment='top', multialignment='left', color='b', fontsize=14)
plt.text(plt.xticks()[0][1], max(plt.yticks()[0])-1, 'Overestimate', rotation=0, horizontalalignment='left', verticalalignment='top', multialignment='left', color='r', fontsize=14)
plt.xlabel('Difference')
plt.savefig('../Paper/Figs/ContainmentHistogram.png')

# Do a true vs containment estimate plot
plt.figure()
f, ax = plt.subplots()
ax.plot([0, 1], [0, 1], ls="--", c=".3")
ax.plot(sorted_true, sorted_containment_estimates)
plt.ylabel('Estimate Jaccard')
plt.xlabel('True Jaccard')
plt.title('Min Hash Via Containment')
axes = plt.gca()
#axes.text(-.2, 1.15, 'b)', horizontalalignment='left', verticalalignment='bottom', fontdict=font)
plt.savefig('../Paper/Figs/ContainmentTrueVsEstimate.png')

# Do a relative error plot
plt.figure()
plt.plot(sorted_true, containment_differences / sorted_true)
axes = plt.gca()
axes.set_ylim([-1, 1])
plt.axhline(0, color='black', linestyle='dashed', linewidth=2)
plt.ylabel('Relative error')
plt.xlabel('True Jaccard index')
plt.savefig('../Paper/Figs/ContainmentRelativeError.png')

# Print out stats and save to file for putting in the paper.
print("Classic Min Hash mean: %f, variance: %f" % (np.mean(differences), np.var(differences)))
print("Containment Min Hash mean: %f, variance: %f" % (np.mean(containment_differences), np.var(containment_differences)))
fid = open(os.path.abspath('../Paper/Data/SyntheticDataClassic.txt'), 'w')
fid.write("$%f\pm%f$" % (np.mean(differences), np.var(differences)))
fid.close()
fid = open(os.path.abspath('../Paper/Data/SyntheticDataContainment.txt'), 'w')
fid.write("$%f\pm%f$" % (np.mean(containment_differences), np.var(containment_differences)))
fid.close()
