#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Created on 28/03/2023
@author: maxime.naour@inrae.fr
"""

import os
import urllib.request
from tqdm import tqdm

genomes = input('Path (if different from the location of the script) and file name containing the genomes to download: ')

# Lire le fichier contenant les informations de téléchargement
with open(genomes, "r") as f:
    lines = f.readlines()[1:] # On ignore la première ligne qui contient les en-têtes

# Pour chaque ligne, télécharger le génome correspondant
for line in tqdm(lines, desc="Downloading genomes", colour="magenta"):
    organism, strain, _, gca, _, gca_url, gcf, _, gcf_url = line.strip().split(", ")

    # Créer un répertoire pour stocker le génome
    directory = f"{organism.replace(' ', '_').replace('.', '')}_{strain.replace(' ', '_').replace('.', '')}"
    os.makedirs(directory, exist_ok=True)

    # Télécharger le génome GCF s'il est disponible, sinon télécharger le génome GCA
    try:
        response = urllib.request.urlopen(gcf_url)
        genome_file = f"{gcf}.fna.gz"
    except:
        response = urllib.request.urlopen(gca_url)
        genome_file = f"{gca}.fna.gz"

    # Écrire le contenu de la réponse dans le fichier génome
    with open(os.path.join(directory, genome_file), "wb") as outfile:
        pbar = tqdm(unit="B", total=int(response.headers["Content-Length"]))
        while True:
            chunk = response.read(1024)
            if not chunk:
                break
            outfile.write(chunk)
            pbar.update(len(chunk))
        pbar.close()
