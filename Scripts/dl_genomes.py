# -*- coding: utf-8 -*-

"""
Created on 28/03/2023
@author: maxime.naour@inrae.fr
"""

import os
import urllib.request
from tqdm import tqdm

genomes = input('Path (if different from the location of the script) and file name containing the genomes to download: ')

# Read the file containing the download information
with open(genomes, "r") as f:
    lines = f.readlines()[1:]  # Ignore the first line containing the headers

# Create the "Genomes" directory to store the downloaded genomes
os.makedirs("Genomes", exist_ok=True)

# Set up the outer progress bar for the overall download process
outer_pbar = tqdm(lines, desc="Downloading genomes", unit=" genome")

# For each line, download the corresponding genome
for line in outer_pbar:
    organism, strain, _, gca, _, gca_url, gcf, _, gcf_url = line.strip().split(", ")

    # Create a directory to store the genome
    directory = f"{organism.replace(' ', '_').replace('.', '')}_{strain.replace(' ', '_').replace('.', '')}"
    os.makedirs(os.path.join("Genomes", directory), exist_ok=True)

    # Download the GCF genome if available, otherwise download the GCA genome
    if gcf != "NA":
        response = urllib.request.urlopen(gcf_url)
        genome_file = f"{gcf}.fna.gz"
    else:
        response = urllib.request.urlopen(gca_url)
        genome_file = f"{gca}.fna.gz"

    # Write the content of the response to the genome file
    with open(os.path.join("Genomes", directory, genome_file), "wb") as outfile:
        while True:
            chunk = response.read(1024)
            if not chunk:
                break
            outfile.write(chunk)

    # Update the outer progress bar
    outer_pbar.update(1)

# Close the outer progress bar
outer_pbar.close()
