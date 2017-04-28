# Since the stupid parallel pool can only be read as a script, here I will read it back in to
# visualize what's going on

import sys
import os, timeit, h5py
import MinHash as MH
import numpy as np
import multiprocessing
from multiprocessing import Pool
from itertools import *
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D

count_range = range(3, 50, 2)

data_dir = '/home/dkoslicki/Dropbox/Repositories/MinHash/data/Mock'
#to_plot = 'Mock_Even_Bacteroides_vulgatus_ATCC_8482.txt'  # This is one in the sample
#to_plot = 'Mock_Even_Bacteroides_vulgatus_CAG_6.txt'  # This is one NOT in the sample

for to_plot in ['Mock_Even_Bacteroides_vulgatus_ATCC_8482_500hash.txt', 'Mock_Even_Bacteroides_vulgatus_CAG_6_500hash.txt']:
	kmer_sizes = []
	hash_sizes = []
	jaccards = []
	jaccard_counts = []
	commons = []
	common_counts = []
	fid = open(os.path.join(data_dir, to_plot), 'r')
	i = 0
	for line in fid.readlines():
		if i == 0:  # Skip the first line
			i += 1
			continue
		line = line.strip()
		(kmer_size, hash_size, jaccard, jaccard_count, common, common_count) = line.split('\t')
		kmer_sizes.append(kmer_size)
		hash_sizes.append(hash_size)
		jaccards.append(jaccard)
		j_c = jaccard_count.replace('(', '')
		j_c = j_c.replace(')', '')
		j_c = np.array(j_c.split(','), dtype=np.float)
		jaccard_counts.append(j_c)
		commons.append(common)
		cc = common_count.replace('(', '')
		cc = cc.replace(')', '')
		cc = np.array(cc.split(','), dtype=np.float)
		common_counts.append(cc)

	fid.close()

	plt.figure()
	x_range = count_range
	line1, = plt.plot(x_range, jaccards, marker='v', label='Jaccards')
	line2, = plt.plot(x_range, [item[0] for item in jaccard_counts], marker='.', label='Jaccard counts 0')
	line3, = plt.plot(x_range, [item[1] for item in jaccard_counts], marker='.', label='Jaccard counts 1')
	line4, = plt.plot(x_range, [item[0]/common_counts[0][0] for item in common_counts], marker='.', label='common count 0 ')
	line5, = plt.plot(x_range, [item[1]/common_counts[0][1] for item in common_counts], marker='.', label='common count 1 ')
	plt.xticks(np.arange(min(x_range), max(x_range)+1, 2.0))
	plt.legend(handler_map={line1: HandlerLine2D(numpoints=4)})
	plt.title(to_plot)

	print(to_plot)
	print('Jaccard AUC: %f' % np.array(jaccards, dtype=np.float).dot(count_range))
	print('Jaccard count 0 AUC: %f' % np.array([item[0] for item in jaccard_counts], dtype=np.float).dot(count_range))
	print('Jaccard count 1 AUC: %f' % np.array([item[1] for item in jaccard_counts], dtype=np.float).dot(count_range))
	print('common counts 0 AUC: %f' % np.array([item[0]/common_counts[0][0] for item in common_counts], dtype=np.float).dot(count_range))
	print('common counts 1 AUC: %f' % np.array([item[1]/common_counts[0][1] for item in common_counts], dtype=np.float).dot(count_range))
