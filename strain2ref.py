#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Created on 28/03/2023
@author: maxime.naour@inrae.fr
"""

"""
This python script allows you to obtain the RefSeq accession number and the "NCBI Taxonomy ID" from a strain name.
This script uses the EDirect tool (language: Perl) to access NCBI (National Center for Biotechnology Information) resources.
"""

import os
import subprocess
from tqdm import tqdm

input_file = input('Path (if different from the location of the script) and file name containing the strain names: ')
output_file = input('Path (if different from the location of the script) and file name to save the corresponding RefSeq accession numbers and taxids: ')

def get_most_recent_refseq_version(refseq_versions):
    refseq_versions_list = refseq_versions.split(';')
    refseq_versions_list = [version.strip() for version in refseq_versions_list]

    most_recent_version = max(refseq_versions_list, key=lambda x: int(x.split('.')[-1]))
    return most_recent_version

def get_refseq_accession_and_taxonomy_id(strain_name):
    esearch = subprocess.Popen(['esearch', '-db', 'biosample', '-query', strain_name], stdout=subprocess.PIPE)
    elink = subprocess.Popen(['elink', '-target', 'nuccore'], stdin=esearch.stdout, stdout=subprocess.PIPE)
    esummary = subprocess.Popen(['esummary'], stdin=elink.stdout, stdout=subprocess.PIPE)
    xtract = subprocess.Popen(['xtract', '-pattern', 'DocumentSummary', '-element', 'AccessionVersion', 'TaxId'], stdin=esummary.stdout, stdout=subprocess.PIPE)
    xtract_output = xtract.communicate()[0].decode().strip().split('\t')

    refseq_accessions = xtract_output[0]
    most_recent_refseq = get_most_recent_refseq_version(refseq_accessions)
    taxonomy_id = xtract_output[-1]

    return most_recent_refseq, taxonomy_id

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    strain_names = infile.read().splitlines()
    for strain_name in tqdm(strain_names, desc="Processing strain names", ncols=100, colour='magenta'):
        refseq_accession, taxonomy_id = get_refseq_accession_and_taxonomy_id(strain_name)
        outfile.write(f"{strain_name}, {refseq_accession}, {taxonomy_id}\n")

header = "Strain Name, RefSeq Accession no, Taxid\n"

with open(output_file, 'r') as outfile:
    content = outfile.read()

with open(output_file, 'w') as outfile:
    outfile.write(header + content)

print('Finished !')
