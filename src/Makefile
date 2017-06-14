dataDir:=../data
paperDir:=../Paper
binDir:=../bin
genSimLoc:=${binDir}/GemSIM_v1.6/GemReads.py
jellyfishLoc:=${binDir}/jellyfish-2.2.3
snapLoc:=${binDir}/snap-aligner

all: ${dataDir}/Viruses/FileNames.txt \
     ${dataDir}/Genomes/FileNames.txt \
     ${dataDir}/SNAP/4539585.3.sorted.r1.fastq \
     ${dataDir}/SNAP/4539585.3.sorted.r2.fastq \
     ${paperDir}/Figs/ContainmentConceptual.png \
     ${paperDir}/Data/SyntheticDataContainment.txt \
     ${dataDir}/Genomes/AllSketches.h5 \
     ${paperDir}/Figs/SimulatedBiologicalData_small.png \
     ${paperDir}/Figs/SimulatedBiologicalData.png \
     ${dataDir}/Viruses/AllSketches.h5 \
     ${paperDir}/Figs/CoveragePlot.png

# Download genSIM
${genSimLoc}:
	mkdir -p ${binDir}
	cd ${binDir} && \
      wget https://sourceforge.net/projects/gemsim/files/GemSIM_v1.6.tar.gz/download && \
      tar xf download && \
	  rm download

# Download jellyfish
${jellyfishLoc}/bin/jellyfish:
	mkdir -p ${binDir}
	cd ${binDir} && \
      wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.3/jellyfish-2.2.3-CentOS6.tar.gz && \
      tar xf jellyfish-2.2.3-CentOS6.tar.gz && \
	  rm jellyfish-2.2.3-CentOS6.tar.gz

# Download jellyfish
${snapLoc}:
	mkdir -p ${binDir}
	cd ${binDir} && \
      wget http://snap.cs.berkeley.edu/downloads/snap-beta.18-linux.tar.gz && \
      tar xf snap-beta.18-linux.tar.gz && \
	  rm snap-beta.18-linux.tar.gz
	chmod +x ${snapLoc}

# Download all of the data I'm using
${dataDir}/SNAP/4539585.3.sorted.r%.fastq:
	wget -nd http://files.cgrb.oregonstate.edu/Koslicki_Lab/MinHashContainment/4539585.3.sorted.r$*.fastq.tar.gz -P ${dataDir}/SNAP/
	tar -xzf ${dataDir}/SNAP/4539585.3.sorted.r$*.fastq.tar.gz -C ${dataDir}/SNAP/
	rm ${dataDir}/SNAP/4539585.3.sorted.r$*.fastq.tar.gz

${dataDir}/Genomes/FileNames.txt:
	wget -nd http://files.cgrb.oregonstate.edu/Koslicki_Lab/MinHashContainment/Genomes.tar.gz -P ${dataDir}
	tar -xzf ${dataDir}/Genomes.tar.gz -C ${dataDir}
	cd ${dataDir}/Genomes && find `pwd` -name "*.fna" > FileNames.txt
	rm ${dataDir}/Genomes.tar.gz

${dataDir}/Viruses/FileNames.txt:
	tar -xzf ${dataDir}/ViralGenomes.tar.gz --strip 1 -C ${dataDir}/Viruses
	cd ${dataDir}/Viruses && find `pwd` -name "*.fna" > FileNames.txt

# Make the conceptual figure
${paperDir}/Figs/ContainmentConceptual.png:
	# Also creates ${paperDir}/Data/ClassicalConceptual.txt
	# Also creates ${paperDir}/Data/ContainmentConceptual.txt
	mkdir -p "${paperDir}/Figs"
	mkdir -p "${paperDir}/Data"
	python ConceptualFigure.py

# Run the synthetic data computations
${paperDir}/Data/SyntheticDataContainment.txt:
	# Also creates:
    #  ${paperDir}/Figs/Differences.png
    #  ${paperDir}/Figs/TrueVsEstimate.png
    #  ${paperDir}/Figs/RelativeError.png
    #  ${paperDir}/Figs/Histogram.png
    #  ${paperDir}/Figs/ContainmentHistogram.png
    #  ${paperDir}/Figs/ContainmentHistogram.png
    #  ${paperDir}/Figs/ContainmentRelativeError.png
    #  ${paperDir}/Data/SyntheticDataClassic.txt
	python JaccardVsContainment.py

# Run the simulated biological data computations
# Creating simulated biological data figures
# Pre-computing bacterial genome sketches (CreateSimulatedMinHashSketches.py)"
${dataDir}/Genomes/AllSketches.h5: ${dataDir}/Genomes/FileNames.txt
	python CreateSimulatedMinHashSketches.py

# Creating small simulated metagenomes and computed jaccard indices and estimates
${paperDir}/Figs/SimulatedBiologicalData_small.png: ${genSimLoc} ${dataDir}/Genomes/FileNames.txt ${dataDir}/Genomes/AllSketches.h5
	# Also creates:
	#   ../data/SimulatedMetagenomes/MinHash_results_small.txt
	#   ../data/SimulatedMetagenomes/ContainmentMinHash_results_small.txt
	#   ../Paper/Data/SimulatedBiologicalData_small_NumGenomes.txt
	#   ../Paper/Data/SimulatedBiologicalData_small_NumReads.txt
	#   ../Paper/Data/SimulatedBiologicalData_small_NumReplicates.txt
	#   ../Paper/Data/SimulatedBiologicalData_small_p.txt
	#   ../Paper/Data/SimulatedBiologicalData_small_ksize.txt
	#   ../Paper/Data/SimulatedBiologicalData_small_rel_size.txt
	python SimulatedBiologicalDataSmall.py ${genSimLoc}

# Creating large simulated metagenomes and computed jaccard indices and estimates
${paperDir}/Figs/SimulatedBiologicalData.png: ${genSimLoc} ${dataDir}/Genomes/FileNames.txt ${dataDir}/Genomes/AllSketches.h5
	# Also creates:
	#   ../data/SimulatedMetagenomes/MinHash_results.txt
	#   ../data/SimulatedMetagenomes/ContainmentMinHash_results.txt
	#   ../Paper/Data/SimulatedBiologicalData_NumGenomes.txt
	#   ../Paper/Data/SimulatedBiologicalData_NumReads.txt
	#   ../Paper/Data/SimulatedBiologicalData_NumReplicates.txt
	#   ../Paper/Data/SimulatedBiologicalData_p.txt
	#   ../Paper/Data/SimulatedBiologicalData_ksize.txt
	#   ../Paper/Data/SimulatedBiologicalData_rel_size.txt
	python SimulatedBiologicalData.py ${genSimLoc}

# Run the real biological data computations
${dataDir}/Viruses/AllSketches.h5: ${dataDir}/Viruses/FileNames.txt
	python CreateVirusesMinHashSketches.py

${paperDir}/Data/MetagenomeTotalKmers.txt: ${jellyfishLoc}/bin/jellyfish ${dataDir}/SNAP/4539585.3.sorted.r1.fastq ${dataDir}/SNAP/4539585.3.sorted.r2.fastq
	# Also creates:
	#   ${dataDir}/MetagenomeBloom.bc
	#   ${dataDir}/MetagenomeBloom.jf
	export LD_LIBRARY_PATH=${jellyfishLoc}/lib:$$LD_LIBRARY_PATH && \
      python MakeMetagenomeBloom.py

${paperDir}/Data/FoundOrganismName.txt: ${paperDir}/Data/MetagenomeTotalKmers.txt ${dataDir}/Viruses/FileNames.txt ${dataDir}/Viruses/AllSketches.h5
	# Also creates:
	#   ../data/MetagenomeViruses.jaccards.txt
	#   ../Paper/Data/FoundOrganismAccession.txt
	#   ../Paper/Data/FoundOrganismFileName.txt
	#   ../Paper/Data/FoundOrganismContainment.txt
	#   ../Paper/Data/FoundOrganismJaccard.txt
	python QueryVirusSketches.py

${paperDir}/Figs/CoveragePlot.png: ${paperDir}/Data/FoundOrganismName.txt ${dataDir}/SNAP/4539585.3.sorted.r1.fastq ${dataDir}/SNAP/4539585.3.sorted.r2.fastq ${snapLoc}
	# Also creates:
	# ${paperDir}/Data/MeanCoverage.txt
	chmod +x MakeCoveragePlot.sh
	./MakeCoveragePlot.sh
