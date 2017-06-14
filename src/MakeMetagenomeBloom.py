import os
import sys
import subprocess
import multiprocessing

num_threads = multiprocessing.cpu_count()
#jellyfish_loc = "/home/dkoslicki/Documents/jellyfish-2.2.3/bin/./jellyfish"
jellyfish_loc = sys.argv[1]
ksize = 21  # k-mer length
p = 0.01
metagenome_file_name_r1 = os.path.abspath("../data/SNAP/4539585.3.sorted.r1.fastq")
metagenome_file_name_r2 = os.path.abspath("../data/SNAP/4539585.3.sorted.r2.fastq")
bloom_out_file = os.path.abspath("../data/MetagenomeBloom.bc")
count_out_file = os.path.abspath("../data/MetagenomeBloom.jf")

# Bloom filter for the actual querying, turns out this actually completely ignores the false positive rate parameter,
# while the jellyfish count appears to respect it, so let's switch to the jellyfish count
#cmd = jellyfish_loc + " bc -s 1G -m " + str(ksize) + " --fpr=" + str(p) + " -C -t " + str(num_threads) + " -o "\
#	+ bloom_out_file + " " + metagenome_file_name_r1 + " " + metagenome_file_name_r2
#test = subprocess.check_output(cmd, shell=True)

# Bloom count filter just for the number of distinct k-mers
cmd = jellyfish_loc + " count -s 100M -m " + str(ksize) + " -C -t " + str(num_threads) + " -o "\
	+ count_out_file + " " + metagenome_file_name_r1 + " " + metagenome_file_name_r2
test = subprocess.check_output(cmd, shell=True)

# Get the number of unique k-mers in the metagenome
cmd = jellyfish_loc + " stats " + count_out_file
res = subprocess.check_output(cmd, shell=True)
num_kmers = int(res.split()[3])
fid = open(os.path.abspath("../Paper/Data/MetagenomeTotalKmers.txt"), 'w')
fid.write("%d" % num_kmers)
fid.close()

# remove count file (only need the bloom filter)
#os.remove(count_out_file)
