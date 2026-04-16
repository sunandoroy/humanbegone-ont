#!/bin/bash

# ======================================================================
# HumanBeGone-ONT Initialization and Setup Script
# Run this once to setup the conda environment and download required databases.
# ======================================================================

# Stop execution if any critical command fails
set -e

# Always execute relative to the script's physical location
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd "$DIR"

echo "=========================================================="
echo "   Initializing HumanBeGone-ONT Environment & Databases"
echo "=========================================================="

# 1. Create the Conda Environment
echo -e "\n[1/6] Setting up the Conda environment from humanbegone-ont.yml..."
if [ -f "humanbegone-ont.yml" ]; then
    conda env create -f humanbegone-ont.yml
    echo "-> Conda environment initialized successfully."
else
    echo "ERROR: humanbegone-ont.yml not found in $DIR"
    exit 1
fi

# Placeholder Zenodo URLs
# Update these links with your actual Zenodo URLs once published.
URL_KRAKEN="https://zenodo.org/records/19591606/files/kraken_index.tar.gz"
URL_REF_FASTA="https://zenodo.org/records/19591606/files/T2T_Unique.fna"
URL_WN_KMER="https://zenodo.org/records/19591606/files/T2T_repetitive_kmers.txt"
URL_TEST_FILES="https://zenodo.org/records/19591606/files/Test-ont.tar.gz"

# 2. Download and Extract Kraken2 Index
echo -e "\n[2/6] Downloading Kraken2 Index..."
wget "$URL_KRAKEN" -O kraken_index.tar.gz
echo "Extracting Kraken2 Index..."
tar -xvzf kraken_index.tar.gz
rm kraken_index.tar.gz

# 3. Download and Extract Reference Fasta
echo -e "\n[3/6] Downloading Reference Fasta..."
wget "$URL_REF_FASTA"

# 4. Download and Extract Winnowmap Kmer Files
echo -e "\n[4/6] Downloading Winnowmap K-mer Table..."
wget "$URL_WN_KMER" 

# 5. Download and Extract Test Files
echo -e "\n[5/6] Downloading Test Files..."
wget "$URL_TEST_FILES" -O test_files.tar.gz
echo "Extracting Test Files..."
tar -xvzf test_files.tar.gz
rm test_files.tar.gz

# 6. Build Minimap2 Index utilizing created Conda Env
echo -e "\n[6/6] Generating Minimap2 Index from the Reference Fasta..."

# Hook into Conda dynamically so bash can natively activate the env
eval "$(conda shell.bash hook)"
conda activate humanbegone-ont

# Modify 'T2T_Unique.fna' below if your downloaded target file is named differently
echo "Building index: minimap2 -d T2T_Human.mmi -k 13 -w 5 T2T_Unique.fna"
minimap2 -d T2T_Human.mmi -k 13 -w 5 T2T_Unique.fna

echo "-> Minimap2 Index Generation Complete."

echo "=========================================================="
echo "   Initialization Complete!"
echo "   Access your databases and test files in:"
echo "   $DIR"
echo "=========================================================="
