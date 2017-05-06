#!/bin/bash
# This script will run all the computations, save the output, automatically populate the LaTeX file, and render the results
# Warning, this can take quite a while to run...

# Need to package and automate the downloading of all the data I'm using

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