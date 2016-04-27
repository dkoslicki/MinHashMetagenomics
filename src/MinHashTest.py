import MinHash as MH

reload(MH)
fid = open('/Users/dkoslicki/Dropbox/Repositories/MinHash/data/test_files.txt', 'r')
file_names = fid.readlines()
fid.close()
file_names = [name.strip() for name in file_names]
CE = MH.CountEstimator(n=5000, ksize=11, input_file_name=file_names[0], save_kmers='y')
CE2 = MH.CountEstimator(n=5000, ksize=11, input_file_name=file_names[1], save_kmers='y')
CE2.jaccard_count(CE)
CE2.jaccard(CE)
CE.jaccard(CE2)

CS = MH.CompositionSketch(n=5000, ksize=11, prefixsize=1, input_file_name=file_names[0])
CS2 = MH.CompositionSketch(n=5000, ksize=11, prefixsize=1, input_file_name=file_names[1])
CS.jaccard_count(CS2)
CS.jaccard(CS2)
CS2.jaccard(CS)


i = 0
for record in screed.open(file_names[0]):
    for kmer in MH.kmers(record.sequence,10):
        if kmer == 'TGGAATTCCA':
            i += 1

print(i)

CE = MH.CountEstimator(n=500, ksize=20, input_file_name=file_names[0], save_kmers='y')
