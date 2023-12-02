#!/save_projet/fine/mnaour/PhAN/Metabarcoding_PhAN/pkg/bin/python3
# -*- coding: utf-8 -*-

"""
Created on 28/02/2023
@author: maxime.naour@inrae.fr

Ce script est utilisé pour extraire des séquences uniques à partir d'un fichier de sortie BLASTN.
Il ajoute également des sauts de ligne entre les séquences qui ne se suivent pas.
"""

import os
import re

# Demande à l'utilisateur le nom du fichier
filename = input("Enter your file's name : ")

lines_list = []
new_data = []

# Vérifie si le fichier existe
if os.path.exists(filename):
    with open(filename, "r") as file:
        lines_list = file.readlines()
else:
    print("The document", filename, "is not detected, make sure you place it in your environment")
    exit()

# Ouvre le fichier de sortie
with open("unique_sequences.txt", 'w') as file_out:
    # Écrit les lignes d'en-tête dans le fichier de sortie
    for line in lines_list:
        if line.startswith("BLASTN") or line.startswith("Reference") or line.startswith("Database") or line.startswith("Query=") or line.startswith("Length="):
            file_out.write(line)

    # Ajoute un saut de ligne après les lignes d'en-tête
    file_out.write("\n")

    # Parcourt toutes les lignes du fichier
    for i in range(len(lines_list)-1):
        lineA = lines_list[i]
        lineB = lines_list[i+1]

        # Si la ligne A contient une requête 'Query_' suivie de chiffres et que la ligne B est un saut de ligne, ajoute la ligne A aux nouvelles données
        if re.search(r'Query_\d+', lineA) and '\n' == lineB:
            new_data.append(lineA)

    # Écrit les nouvelles données dans le fichier de sortie
    for line in new_data:
        file_out.write(line)

# Section pour ajouter des sauts de ligne entre les séquences qui ne se suivent pas
with open("unique_sequences.txt", 'r') as file:
    lines_list = file.readlines()

new_data = []
start_adding_newlines = False

# Parcourt toutes les lignes du fichier
for i in range(0, len(lines_list)-1):
    lineA = lines_list[i]
    lineB = lines_list[i+1]

    # Utilise une expression régulière pour identifier toutes les occurrences de 'Query_' suivies de chiffres
    if re.search(r'Query_\d+', lineA):
        start_adding_newlines = True

    if start_adding_newlines:
        if len(lineA.split()) > 3 and len(lineB.split()) > 1:
            end_num_A = int(lineA.split()[3])
            start_num_B = int(lineB.split()[1])

            if end_num_A + 1 != start_num_B:
                new_data.append(lineA)
                new_data.append('\n')
            else:
                new_data.append(lineA)

# Écrit les nouvelles données dans le fichier de sortie
with open("unique_sequences.txt", 'w') as file_out:
    for line in new_data:
        file_out.write(line)

print("Done. The unique sequences have been extracted to unique_sequences.txt.")
