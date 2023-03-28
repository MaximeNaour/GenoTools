#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Created on 28/03/2023
@author: maxime.naour@inrae.fr
"""

"""
This python script allows you to obtain the full name of the strain, the RefSeq accession number, and the "NCBI taxonomy ID" from a simple Genbank accession number.
This script uses the EDirect tool (language: Perl) to access NCBI (National Center for Biotechnology Information) resources.
"""

import os
import subprocess
import re
from tqdm import tqdm

input_file = input('Path (if different from the location of the script) and file name containing the GenBank accessions: ')
output_file = input('Path (if different from the location of the script) and file name to save the corresponding organism and strain names: ')

# Grouper les queries en lots de 10 pour gagner du temps
batch_size = 10

def get_refseq_accession_and_taxonomy_id(organism, strain):
    query = f'"{organism}"[Organism] AND "{strain}"[All Fields]'
    esearch = subprocess.Popen(['esearch', '-db', 'assembly', '-query', query], stdout=subprocess.PIPE)
    elink = subprocess.Popen(['elink', '-target', 'nuccore', '-name', 'assembly_nuccore_refseq'], stdin=esearch.stdout, stdout=subprocess.PIPE)
    esummary = subprocess.Popen(['esummary'], stdin=elink.stdout, stdout=subprocess.PIPE)
    xtract = subprocess.Popen(['xtract', '-pattern', 'DocumentSummary', '-element', 'AccessionVersion', 'TaxId'], stdin=esummary.stdout, stdout=subprocess.PIPE)
    xtract_output = xtract.communicate()[0].decode().strip().split('\t')
    
    refseq_accession = xtract_output[0]
    taxonomy_id = xtract_output[-1]

    return refseq_accession, taxonomy_id


with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    accessions = infile.read().splitlines()
    for acc in tqdm(accessions, desc="Processing accessions", ncols=100, colour='magenta'):
        esearch = subprocess.Popen(['esearch', '-db', 'nuccore', '-query', acc], stdout=subprocess.PIPE)
        esummary = subprocess.Popen(['esummary'], stdin=esearch.stdout, stdout=subprocess.PIPE)
        xtract = subprocess.Popen(['xtract', '-pattern', 'DocumentSummary', '-element', 'Title'], stdin=esummary.stdout, stdout=subprocess.PIPE)
        organism_name = xtract.communicate()[0].decode().strip()

        strain_pattern = re.compile(r'strain\s+([^\s,;]+)')
        strain_match = strain_pattern.search(organism_name)

        if strain_match:
            strain = strain_match.group(1)
            organism = organism_name.split(' strain ')[0]
        else:
            organism = organism_name
            strain = ''

        refseq_accession, taxonomy_id = get_refseq_accession_and_taxonomy_id(organism, strain)

        outfile.write(f"{acc}, {organism.split(',')[0]}, {strain}, {refseq_accession}, {taxonomy_id}\n")

# Supprimer les lignes exédentaires
with open(output_file, 'r') as outfile:
    lines = outfile.read().splitlines()

modified_lines = []
for line in lines:
    # Vérifier si le premier mot de la ligne se trouve dans la liste des accessions
    first_word = line.split(',')[0]
    if first_word in accessions:
        modified_lines.append(line.rstrip(',') + '\n')

with open(output_file, 'w') as outfile:
    outfile.writelines(modified_lines)
    
header = "Accession no, Organism, Strain, RefSeq no, taxid\n"

with open(output_file, 'r') as outfile:
    content = outfile.read()

with open(output_file, 'w') as outfile:
    outfile.write(header + content)

print('Finished !')
