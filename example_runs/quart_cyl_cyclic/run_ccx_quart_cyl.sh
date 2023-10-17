#!/bin/bash

echo " "
echo "Initializing Anaconda base environment..."

## ---------- INITIALIZE ANACONDA ENVIRONMENT FOR LATER ---------- 

## My path to the Anaconda initialization script is in "/opt/"
## Another common installation directory is, "~/anaconda3/"
## Include the trailing slash in this filepath definition.
CONDA_PATH="/opt/anaconda3/" 

## Filepath to the key executable within the Anaconda3 directory. This should not
## require any changes with a default installation.
CONDA_SH_SUBPATH="etc/profile.d/conda.sh"

## Now actually execute the Anaconda intialization script
source "$CONDA_PATH$CONDA_SH_SUBPATH"

## And activate the "base" environment. I am assuming that the "base" Anaconda
## environment contains python3, numpy, scipy, and pandas.
conda activate base

echo "Success!"
echo " "


## ---------- BEGIN CALCULIX SIMULATION ---------- 

echo "Starting CalculiX simulation..."

# Set this to the number of threads to run in parallel for the CCX simulation
export OMP_NUM_THREADS=4

# Filepath to the CalculiX executable you wish to use
CCX_BIN="../bin/ccx_2.20_MT_SLS_umat"

## Filepath of the simulation input file without the extension. If you created
## an input file called, "my_job.inp", then set this variable equal to "my_job".
## In general, it will be safest to execute this script in the same directory
## as the input file itself.
INPUT_FILE="./quart_cyl_cycle_compress_umat"

## Run the CCX simulation
$CCX_BIN $INPUT_FILE >> stdout.log

echo "Success!"
echo " "


## ---------- BEGIN PYTHON POST-PROCESSING ---------- 

echo "Starting CalculiX post-processing..."
echo " "

## Run the CCX utility script that parses the .dat file. This will write some
## text files to the hard drive which are slightly easier to parse than the
## basic .dat files.
python3 ../../src/dat2txt.py $INPUT_FILE

## Run the post-analysis script
python3 ../../src/analyze_viscoelast_output.py
echo " "