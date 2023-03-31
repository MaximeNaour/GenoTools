#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Created on 31/03/2023
@author: maxime.naour@inrae.fr
"""

# pip install biopython

import csv
from Bio import Entrez

Entrez.email = "maxime.naour@inrae.fr"

def search_ncbi(species, strain_info):
    try:
        query = f"{species} {strain_info}"
        handle = Entrez.esearch(db="taxonomy", term=query)
        record = Entrez.read(handle)
        handle.close()

        if record["Count"] == "1":
            tax_id = record["IdList"][0]
            handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
            records = Entrez.read(handle)
            handle.close()

            scientific_name = records[0]["ScientificName"]
            if strain_info not in scientific_name:
                return f"{scientific_name} {strain_info}"
            else:
                return scientific_name
        else:
            return None
    except Exception as e:
        print(f"Error: {e}")
        return None

def update_bacteria_names(input_file, output_file):
    with open(input_file, "r") as input_csv, open(output_file, "w", newline='') as output_csv:
        reader = csv.reader(input_csv, delimiter=' ')
        
        writer = csv.writer(output_csv, delimiter=',')
        writer.writerow(['Species', 'Strain Info', 'NCBI Name'])

        for row in reader:
            species = row[0] + ' ' + row[1]
            strain_info = ' '.join(row[2:])

            ncbi_name = search_ncbi(species, strain_info)

            if ncbi_name:
                writer.writerow([species, strain_info, ncbi_name])
            else:
                writer.writerow([species, strain_info, 'Not found'])

input_file = input('Path (if different from the location of the script) and file name containing the organism and strain names: ')
output_file = input('Path (if different from the location of the script) and file name to save the correction of organism and strain names: ')

update_bacteria_names(input_file, output_file)
