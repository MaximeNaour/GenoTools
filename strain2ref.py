#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Created on 28/03/2023
@author: maxime.naour@inrae.fr
"""

"""
This python script allows you to obtain the most recent GenBank assembly accession number (GCA accession number) URL and RefSeq assembly accession number (GCF accession number) URL from an organism and strain name.
This script uses the EDirect tool (language: Perl) to access NCBI (National Center for Biotechnology Information) resources.
"""

import os
import re
import subprocess
from tqdm import tqdm

input_file = input('Path (if different from the location of the script) and file name containing the organism and strain names: ')
output_file = input('Path (if different from the location of the script) and file name to save the corresponding GenBank assembly accession URLs and RefSeq assembly accession URLs: ')

header_input = input('Is there a header in the input file? (yes/no): ').lower()

### Unused in this script !
def get_most_recent_refseq_version(refseq_versions):
    refseq_versions_list = refseq_versions.split(';')
    refseq_versions_list = [version.strip() for version in refseq_versions_list]

    most_recent_version = max(refseq_versions_list, key=lambda x: int(x.split('.')[-1]))
    return most_recent_version
###

def get_gca_url_and_taxonomy(organism, strain):
    query = f'"{organism}"[Organism] AND "{strain}"[All Fields] AND latest[filter]'
    esearch = subprocess.Popen(['esearch', '-db', 'assembly', '-query', query], stdout=subprocess.PIPE)
    esummary = subprocess.Popen(['esummary'], stdin=esearch.stdout, stdout=subprocess.PIPE)
    xtract = subprocess.Popen(['xtract', '-pattern', 'DocumentSummary', '-element', 'FtpPath_GenBank'], stdin=esummary.stdout, stdout=subprocess.PIPE)
    ftp_paths = xtract.communicate()[0].decode().strip().splitlines()

    if not ftp_paths:
        print(f"WARNING: No GenBank assembly accession URL found for {organism}, {strain}")
        return 'NA', 'NA', 'NA', 'NA'

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
    gca_assembly_version = gca_id_and_version.replace(gca_id + '_', '')
    
    esearch_tax = subprocess.Popen(['esearch', '-db', 'taxonomy', '-query', organism], stdout=subprocess.PIPE)
    esummary_tax = subprocess.Popen(['esummary'], stdin=esearch_tax.stdout, stdout=subprocess.PIPE)
    xtract_tax = subprocess.Popen(['xtract', '-pattern', 'DocumentSummary', '-element', 'TaxId'], stdin=esummary_tax.stdout, stdout=subprocess.PIPE)
    taxonomy_id = xtract_tax.communicate()[0].decode().strip()
    
    return gca_id, gca_assembly_version, latest_gca_url, taxonomy_id

def get_gcf_url(organism, strain):
    query = f'"{organism}"[Organism] AND "{strain}"[All Fields] AND latest[filter]'
    esearch = subprocess.Popen(['esearch', '-db', 'assembly', '-query', query], stdout=subprocess.PIPE)
    esummary = subprocess.Popen(['esummary'], stdin=esearch.stdout, stdout=subprocess.PIPE)
    xtract = subprocess.Popen(['xtract', '-pattern', 'DocumentSummary', '-element', 'FtpPath_RefSeq'], stdin=esummary.stdout, stdout=subprocess.PIPE)
    ftp_paths = xtract.communicate()[0].decode().strip().splitlines()

    if not ftp_paths:
        print(f"WARNING: No RefSeq assembly accession URL found for {organism}, {strain}")
        return 'NA', 'NA', 'NA'

    gcf_urls = []
    for ftp_path in ftp_paths:
        gcf_file_name = ftp_path.split('/')[-1] + '_genomic.fna.gz'
        gcf_url = f"{ftp_path}/{gcf_file_name}"
        gcf_urls.append(gcf_url)

    # Sort GCF URLs by assembly version number (highest first)
    gcf_urls.sort(key=lambda x: int(re.search(r'GCF_\d+\.(\d+)', x).group(1)), reverse=True)

    latest_gcf_url = gcf_urls[0]
    gcf_id_and_version = latest_gcf_url.split('/')[-2]
    gcf_id = re.match(r'(GCF_\d+\.\d+)', gcf_id_and_version).group(1)
    gcf_assembly_version = gcf_id_and_version.replace(gcf_id + '_', '')
    
    return gcf_id, gcf_assembly_version, latest_gcf_url

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    organism_strain_pairs = infile.read().splitlines()

    # Add these lines to skip the header if the user answers "yes"
    if header_input == 'yes':
        organism_strain_pairs.pop(0)

    for organism_strain in tqdm(organism_strain_pairs, desc="Processing organism and strain names", ncols=100, colour='magenta'):
        # Utilisez une expression régulière pour diviser la chaîne avec ou sans espace après la virgule
        organism_strain_split = re.split(r',\s?', organism_strain)
    
        # Vérifiez si la chaîne a été correctement divisée en deux éléments
        if len(organism_strain_split) != 2:
            print(f"WARNING: Invalid format for line '{organism_strain}', skipping.")
            continue
    
        organism, strain = organism_strain_split
        gca_id, gca_assembly_version, gca_url, taxonomy_id = get_gca_url_and_taxonomy(organism, strain)
        gcf_id, gcf_assembly_version, gcf_url = get_gcf_url(organism, strain)
        outfile.write(f"{organism}, {strain}, {taxonomy_id}, {gca_id}, {gca_assembly_version}, {gca_url}, {gcf_id}, {gcf_assembly_version}, {gcf_url}\n")

header = "Organism, Strain, NCBI Taxid, GenBank Assembly Accession, Assembly Version, GenBank Assembly Accession URL, RefSeq Assembly Accession, RefSeq Version, RefSeq Assembly Accession URL\n"

with open(output_file, 'r') as outfile:
    content = outfile.read()

with open(output_file, 'w') as outfile:
    outfile.write(header + content)

print('Finished !')
