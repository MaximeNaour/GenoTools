#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Created on 19/04/2023
@author: maxime.naour@inrae.fr
"""

import csv
import os
import subprocess
from pathlib import Path

def read_input_file(input_file):
    strains_data = []
    with open(input_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            strains_data.append(row)
    return strains_data

def edirect_query(refseq_accession):
    query = f'esearch -db nuccore -query "{refseq_accession}" | efetch -format gb'
    output = subprocess.check_output(query, shell=True, text=True)
    return output

def parse_genes_functions(gb_output):
    genes_functions = []
    lines = gb_output.split('\n')
    for line in lines:
        if line.startswith("     gene            "):
            gene = line.split()[-1]
        if line.startswith("                     /locus_tag="):
            locus_tag = line.split('=')[-1].strip('"')
        if line.startswith("                     /product="):
            product = line.split('=')[-1].strip('"')
            genes_functions.append((gene, locus_tag, product))
    return genes_functions

def write_output_file(strain, genes_functions):
    output_file = f"{strain}_genes_functions.csv"
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['gene', 'locus_tag', 'function']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for gene, locus_tag, function in genes_functions:
            writer.writerow({'gene': gene, 'locus_tag': locus_tag, 'function': function})

def main():
    input_file = "input.csv"
    strains_data = read_input_file(input_file)

    for strain_data in strains_data:
        refseq_accession = strain_data["RefSeq Assembly Accession"]
        strain = strain_data["Strain"]

        print(f"Processing strain {strain} with RefSeq Accession {refseq_accession}")

        gb_output = edirect_query(refseq_accession)
        genes_functions = parse_genes_functions(gb_output)
        write_output_file(strain, genes_functions)

        print(f"Completed strain {strain}")

if __name__ == "__main__":
    main()
