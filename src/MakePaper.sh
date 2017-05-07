#!/bin/bash
# This script will run all the computations, save the output, automatically populate the LaTeX file, and render the results
# Warning, this can take quite a while to run...
dataDir="../data"
mkdir -p ${dataDir}/SNAP
mkdir -p ${dataDir}/Viruses
mkdir -p "../Paper"
mkdir -p "../Paper/Figs"
mkdir -p "../Paper/Data"

# Need to package and automate the downloading of all the data I'm using
echo "Downloading Data"
wget -nd http://files.cgrb.oregonstate.edu/Koslicki_Lab/MinHashContainment/4539585.3.sorted.r1.fastq.tar.gz -P ${dataDir}/SNAP/
wget -nd http://files.cgrb.oregonstate.edu/Koslicki_Lab/MinHashContainment/4539585.3.sorted.r2.fastq.tar.gz -P ${dataDir}/SNAP/
wget -nd http://files.cgrb.oregonstate.edu/Koslicki_Lab/MinHashContainment/ViralGenomes.tar.gz -P ${dataDir}
wget -nd http://files.cgrb.oregonstate.edu/Koslicki_Lab/MinHashContainment/Genomes.tar.gz -P ${dataDir}
echo "Decompressing downloaded data"
tar -xzf ${dataDir}/SNAP/4539585.3.sorted.r1.fastq.tar.gz -C ${dataDir}/SNAP/
tar -xzf ${dataDir}/SNAP/4539585.3.sorted.r2.fastq.tar.gz -C ${dataDir}/SNAP/
tar -xzf ${dataDir}/ViralGenomes.tar.gz -C ${dataDir}
mv ${dataDir}/Genomes/* ${dataDir}/Viruses
rmdir ${dataDir}/Genomes
tar -xzf ${dataDir}/Genomes.tar.gz -C ${dataDir}

# Make the conceptual figure
echo "Creating conceptual figure"
python ConceptualFigure.py

# Run the synthetic data computations
echo "Creating synthetic data figures"
python JaccardVsContainment.py

# Run the simulated biological data computations
echo "Creating simulated biological data figures"
python CreateSimulatedMinHashSketches.py
python SimulatedBiologicalDataSmall.py
python SimulatedBiologicalData.py

# Run the real biological data computations
echo "Creating real data figure(s)"
python CreateVirusesMinHashSketches.py
# NEED TO ADD SCRIPT TO COMPUTE THE JELLYFISH BLOOM FILTER metagenome_bloom_filter
python QueryVirusSketches.py  # Need to make it automatically spit out the top guy
# Need to automate the downloading of this top guy and the poplation of the MakeCoveragePlot.sh target
chmod +x MakeCoveragePlot.sh
./MakeCoveragePlot.sh

# Compile the LaTeX paper
pdflatex ../Paper/ImprovedMinHashForMetagenomics.tex
bibtex ../Paper/ImprovedMinHashForMetagenomics.aux
pdflatex ../Paper/ImprovedMinHashForMetagenomics.tex
pdflatex ../Paper/ImprovedMinHashForMetagenomics.tex