#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Created on 28/03/2023
@author: maxime.naour@inrae.fr
"""

"""
This python script allows you to obtain the most recent GenBank assembly accession number (GCA accession number) URL from an organism and strain name.
This script uses the EDirect tool (language: Perl) to access NCBI (National Center for Biotechnology Information) resources.
"""

import os
import re
import subprocess
from tqdm import tqdm

input_file = input('Path (if different from the location of the script) and file name containing the organism and strain names: ')
output_file = input('Path (if different from the location of the script) and file name to save the corresponding GenBank assembly accession URLs: ')

### Unused in this script !
def get_most_recent_refseq_version(refseq_versions):
    refseq_versions_list = refseq_versions.split(';')
    refseq_versions_list = [version.strip() for version in refseq_versions_list]

    most_recent_version = max(refseq_versions_list, key=lambda x: int(x.split('.')[-1]))
    return most_recent_version
###

def get_gca_url(organism, strain):
    query = f'"{organism}"[Organism] AND "{strain}"[All Fields] AND latest[filter]'
    esearch = subprocess.Popen(['esearch', '-db', 'assembly', '-query', query], stdout=subprocess.PIPE)
    esummary = subprocess.Popen(['esummary'], stdin=esearch.stdout, stdout=subprocess.PIPE)
    xtract = subprocess.Popen(['xtract', '-pattern', 'DocumentSummary', '-element', 'FtpPath_GenBank'], stdin=esummary.stdout, stdout=subprocess.PIPE)
    ftp_paths = xtract.communicate()[0].decode().strip().splitlines()

    if not ftp_paths:
        print(f"WARNING: No GenBank assembly accession URL found for {organism}, {strain}")
        return '', '', ''

    gca_urls = []
    for ftp_path in ftp_paths:
        gca_file_name = ftp_path.split('/')[-1] + '_genomic.fna.gz'
        gca_url = f"{ftp_path}/{gca_file_name}"
        gca_urls.append(gca_url)

    # Sort GCA URLs by assembly version number (highest first)
    gca_urls.sort(key=lambda x: int(re.search(r'GCA_\d+\.(\d+)', x).group(1)), reverse=True)

    latest_gca_url = gca_urls[0]
    gca_id_and_version = latest_gca_url.split('/')[-2]
    gca_id = re.match(r'(GCA_\d+\.\d+)', gca_id_and_version).group(1)
    assembly_version = gca_id_and_version.replace(gca_id + '_', '')
    
    return gca_id, assembly_version, latest_gca_url

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    organism_strain_pairs = infile.read().splitlines()
    for organism_strain in tqdm(organism_strain_pairs, desc="Processing organism and strain names", ncols=100, colour='magenta'):
        organism, strain = organism_strain.split(', ')
        gca_id, assembly_version, gca_url = get_gca_url(organism, strain)
        outfile.write(f"{organism}, {strain}, {gca_id}, {assembly_version}, {gca_url}\n")

header = "Organism, Strain, GenBank Assembly Accession, Assembly Version, GenBank Assembly Accession URL\n"

with open(output_file, 'r') as outfile:
    content = outfile.read()

with open(output_file, 'w') as outfile:
    outfile.write(header + content)

print('Finished !')
