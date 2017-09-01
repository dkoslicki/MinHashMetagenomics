#!/bin/bash
# This script will run all the computations, save the output, automatically populate the LaTeX file, and render the results
# Warning, this can take quite a while to run...
srcDir=`pwd`
dataDir="../data"
genSimLoc="/home/dkoslicki/Documents/GemSIM_v1.6/GemReads.py"
python2Loc="/usr/bin/python"  # hash function not stable between python versions, and stupid byte encoding issues with python 3
mkdir -p ${dataDir}/SNAP
mkdir -p ${dataDir}/Viruses
mkdir -p "../Paper"
mkdir -p "../Paper/Figs"
mkdir -p "../Paper/Data"
mkdir -p ${dataDir}/SimulatedMetagenomes

# Download all of the data I'm using
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
# Create the file names file
cd ${dataDir}/Genomes
find `pwd` -name "*.fna" > FileNames.txt
cd ${srcDir}
cd ${dataDir}/Viruses
find `pwd` -name "*.fna" > FileNames.txt
cd ${srcDir}

# Make the conceptual figure
echo "Creating conceptual figure (ConceptualFigure.py)"
python ConceptualFigure.py

# Make the figures of the probability of deviation and number of hash functions
python DeltaHashesConfPlotScript.py

# Make the figure of number of hash functions as a function of |B|/|A|
python JaccardIndexNumOfHashes.py

# Run the synthetic data computations
echo "Creating synthetic data figures (JaccardVsContainment.py)"
python JaccardVsContainment.py

# Run the simulated biological data computations
echo "Creating simulated biological data figures"
echo "Pre-computing bacterial genome sketches (CreateSimulatedMinHashSketches.py)"
${python2Loc} CreateSimulatedMinHashSketches.py
echo "Creating small simulated metagenomes and computed jaccard indicies and estimates (SimulatedBiologicalDataSmall.py)"
${python2Loc} SimulatedBiologicalDataSmall.py ${genSimLoc} ${python2Loc}
echo "Creating large simulated metagenomes and computed jaccard indicies and estimates (SimulatedBiologicalData.py)"
${python2Loc} SimulatedBiologicalData.py ${genSimLoc} ${python2Loc}

# Run the real biological data computations
echo "Creating real data figure(s)"
python CreateVirusesMinHashSketches.py
python MakeMetagenomeBloom.py
python QueryVirusSketches.py
chmod +x MakeCoveragePlot.sh
./MakeCoveragePlot.sh

# Trim whitespace from all figures
ls ../Paper/Figs/*.png | xargs -I{} convert {} -trim {}

# Compile the LaTeX paper
pdflatex ../Paper/ImprovedMinHashForMetagenomics.tex
bibtex ../Paper/ImprovedMinHashForMetagenomics.aux
pdflatex ../Paper/ImprovedMinHashForMetagenomics.tex
pdflatex ../Paper/ImprovedMinHashForMetagenomics.tex