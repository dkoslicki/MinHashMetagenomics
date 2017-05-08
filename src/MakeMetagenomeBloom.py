import os
import sys
import subprocess
import multiprocessing

num_threads = multiprocessing.cpu_count()
jellyfish_loc = "/home/dkoslicki/Documents/jellyfish-2.2.3/bin/./jellyfish"
ksize = 21  # k-mer length
p = 0.01
metagenome_file_name_r1 = os.path.abspath("../data/SNAP/4539585.3.sorted.r1.fastq")
metagenome_file_name_r2 = os.path.abspath("../data/SNAP/4539585.3.sorted.r2.fastq")
bloom_out_file = os.path.abspath("../data/MetagenomeBloom.bc")

cmd = jellyfish_loc + " bc -s 1G -m " + str(ksize) + " --fpr=" + str(p) + " -C -t " + str(num_threads) + " -o "\
	+ bloom_out_file + " " + metagenome_file_name_r1 + " " + metagenome_file_name_r2

test = subprocess.check_output(cmd, shell=True)

