#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Created on 28/03/2023
@author: maxime.naour@inrae.fr
"""

import os
import sys
import glob
import gzip
import tempfile
from Bio import SeqIO
from subprocess import call

def decompress_gz_file(input_file):
    with gzip.open(input_file, 'rb') as f_in:
        with tempfile.NamedTemporaryFile(mode="wb", suffix=".fna", delete=False) as f_out:
            f_out.writelines(f_in)
            return f_out.name

def annotate_genomes(input_files):
    for input_file in input_files:
        decompressed_file = decompress_gz_file(input_file)
        base_name = os.path.splitext(os.path.splitext(os.path.basename(decompressed_file))[0])[0]
        input_dir = os.path.dirname(input_file)
        output_dir = os.path.join(input_dir, base_name)

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        call(["/usr/local/genome/Anaconda3/envs/prokka-1.14.6/bin/prokka", decompressed_file, "--outdir", output_dir, "--prefix", base_name, "--force"])

        os.remove(decompressed_file)

def parse_prokka_gff(gff_file):
    with open(gff_file, "r") as file:
        for line in file:
            if line.startswith("##sequence-region"):
                continue

            columns = line.strip().split("\t")

            if len(columns) != 9:
                continue

            seq_id, source, feature_type, start, end, score, strand, phase, attributes = columns
            attributes_dict = {key: value for key, value in [attribute.split("=") for attribute in attributes.split(";")]}

            if feature_type == "CDS":
                gene_id = attributes_dict["ID"]
                product = attributes_dict["product"]
                print(f"{gene_id}\t{product}")

def main():
    input_folder = "Genomes"
    input_files = []

    # Recherche des fichiers ".fna.gz" dans le répertoire "Genomes" et ses sous-dossiers
    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if file.endswith(".fna.gz"):
                input_files.append(os.path.join(root, file))

    print(f"Compressed fasta files found : {len(input_files)}")
    for file in input_files:
        print(file)

    annotate_genomes(input_files)

    for input_file in input_files:
        decompressed_file = input_file[:-3]
        base_name = os.path.splitext(os.path.splitext(os.path.basename(decompressed_file))[0])[0]
        input_dir = os.path.dirname(input_file)
        gff_file = os.path.join(input_dir, base_name, f"{base_name}.gff")
        
        print(f"\nGenes in {input_file}:")
        parse_prokka_gff(gff_file)

if __name__ == "__main__":
    main()
