#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Created on 28/03/2023
@author: maxime.naour@inrae.fr
"""

"""
This python script allows you to obtain the full name of the strain, the latest RefSeq GCF accession number, the RefSeq accession number, and the "NCBI taxo
nomy ID" from a simple Genbank accession number.
This script uses the EDirect tool (language: Perl) to access NCBI (National Center for Biotechnology Information) resources.
"""

import os
import requests
import csv
import re
from io import StringIO
from pathlib import Path
from Bio import Entrez, SeqIO

###### User Configuration ######
input_file = "input.csv"
output_file = "compare_outputs.txt"
email = "maxime.naour@inrae.fr"
############# END ##############

# Configurer l'e-mail pour Entrez
Entrez.email = email

# Supprime le fichier de sortie s'il existe déjà
if os.path.isfile(output_file):
    os.remove(output_file)

# Convertir les fins de ligne du fichier d'entrée en format Unix
with open(input_file, 'r') as f:
    content = f.read()
    content = content.replace('\r\n', '\n')
with open(input_file, 'w') as f:
    f.write(content)

# Parcourir chaque ligne du fichier d'entrée
with open(input_file, 'r') as f:
    input_data = csv.reader(f)
    next(input_data)  # Ignorer l'en-tête

    count = 0
    for refseq_assembly, ncbi_reference_sequence in input_data:
        count += 1
        print(f"Traitement de la ligne {count}: {refseq_assembly}, {ncbi_reference_sequence}")

        # Récupérer l'URL du fichier assembly_summary.txt pour l'organisme d'intérêt
        response = requests.get("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt")
        assembly_summary = response.text
        summary_url = ""
        for line in assembly_summary.splitlines():
            fields = line.split("\t")
            if fields[0] == refseq_assembly:
                summary_url = fields[19]
                break

        # Récupérer le nom du fichier FASTA à partir de l'URL du résumé
        fasta_file_name = Path(summary_url + "_genomic.fna.gz").name

        # Télécharger le fichier FASTA
        fasta_url = f"{summary_url}/{fasta_file_name}"
        response = requests.get(fasta_url)
        with open(fasta_file_name, 'wb') as f:
            f.write(response.content)

        # Décompresser le fichier FASTA
        os.system(f"gzip -d {fasta_file_name}")

        # Lire le fichier FASTA décompressé
        fasta_file_content = ""
        with open(fasta_file_name[:-3], 'r') as f:
            fasta_file_content = f.read()

        # Récupérer les séquences pour les deux identifiants de séquence
        sequence1 = fasta_file_content
        handle = Entrez.efetch(db="nuccore", id=ncbi_reference_sequence, rettype="fasta", retmode="text")
        sequence2 = handle.read()

        # Supprimer les informations d'en-tête des séquences
        sequence1_clean = re.sub(r">.*\n", "", sequence1).replace('\n', '').replace('\r', '')
        sequence2_clean = re.sub(r">.*\n", "", sequence2).replace('\n', '').replace('\r', '')

        # Compter le nombre de nucléotides pour chaque séquence
        sequence1_length = len(sequence1_clean)
        sequence2_length = len(sequence2_clean)

        # Comparer les nombres de nucléotides et écrire le résultat dans le fichier de sortie
        with open(output_file, 'a') as f:
            f.write(f"Comparaison pour {refseq_assembly} et {ncbi_reference_sequence}\n")
            if sequence1_length == sequence2_length:
                f.write(f"Le nombre total de nucléotides est identique pour {refseq_assembly} et {ncbi_reference_sequence}.\n")
                f.write(f"{refseq_assembly} : {sequence1_length} nucléotides\n")
                f.write(f"{ncbi_reference_sequence} : {sequence2_length} nucléotides\n")
            else:
                f.write(f"Le nombre total de nucléotides n'est pas identique pour {refseq_assembly} et {ncbi_reference_sequence}.\n")
                f.write(f"{refseq_assembly} : {sequence1_length} nucléotides\n")
                f.write(f"{ncbi_reference_sequence} : {sequence2_length} nucléotides\n")

            # Ajouter un saut de ligne dans le fichier de sortie pour séparer chaque comparaison
            f.write("\n")

        # Supprimer le fichier FASTA téléchargé
        os.remove(fasta_file_name[:-3])
