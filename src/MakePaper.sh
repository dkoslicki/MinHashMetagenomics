#!/bin/bash
# This script will run all the computations, save the output, automatically populate the LaTeX file, and render the results
# Warning, this can take quite a while to run...
dataDir="../data"
mkdir -p ${dataDir}/SNAP
mkdir -p ${dataDir}/Viruses

# Need to package and automate the downloading of all the data I'm using
wget -nd http://files.cgrb.oregonstate.edu/Koslicki_Lab/MinHashContainment/4539585.3.sorted.r1.fastq.tar.gz -P ${dataDir}/Snap/
wget -nd http://files.cgrb.oregonstate.edu/Koslicki_Lab/MinHashContainment/4539585.3.sorted.r2.fastq.tar.gz -P ${dataDir}/Snap/
wget -nd http://files.cgrb.oregonstate.edu/Koslicki_Lab/MinHashContainment/ViralGenomes.tar.gz -P ${dataDir}
tar -xzf ${dataDir}/Snap/4539585.3.sorted.r1.fastq.tar.gz
tar -xzf ${dataDir}/Snap/4539585.3.sorted.r2.fastq.tar.gz
tar -xzf ${dataDir}/ViralGenomes.tar.gz
mv ${dataDir}/Genomes/* ${dataDir}/Viruses
rmdir ${dataDir}/Genomes

# Make the conceptual figure
python ConceptualFigure.py

# Run the synthetic data computations
python JaccardVsContainment.py

# Run the simulated biological data computations
python CreateSimulatedMinHashSketches.py
python SimulatedBiologicalDataSmall.py
python SimulatedBiologicalData.py

# Run the real biological data computations
python CreateVirusesMinHashSketches.py
# NEED TO ADD SCRIPT TO COMPUTE THE JELLYFISH BLOOM FILTER metagenome_bloom_filter
python QueryVirusSketches.py  # Need to make it automatically spit out the top guy
# Need to automate the downloading of this top guy and the poplation of the MakeCoveragePlot.sh target
chmod +x MakeCoveragePlot.sh
./MakeCoveragePlot.sh