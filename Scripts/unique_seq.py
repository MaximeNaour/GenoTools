#!/save_projet/fine/mnaour/PhAN/Metabarcoding_PhAN/pkg/bin/python3
# -*- coding: utf-8 -*-

"""
Created on 28/02/2023
@author: maxime.naour@inrae.fr

Ce script est utilisé pour extraire des séquences uniques à partir d'un fichier de sortie BLASTN.
Il ajoute également des sauts de ligne entre les séquences qui ne se suivent pas.
"""

import os

# Demande à l'utilisateur le nom du fichier et le nombre de requêtes
filename = input("Enter your file's name : ")
num = input("Enter the number of queries : ")

lines_list = []
new_data = []

# Vérifie si le fichier existe
if os.path.exists(filename):

    # Ouvre le fichier et lit toutes les lignes
    with open(filename,"r") as file :
        lines_list = file.readlines()

else:
    # Si le fichier n'existe pas, affiche un message d'erreur
    print("The document", filename, "is not detected, make sure you place it in your environment")

# Ouvre le fichier de sortie
with open("unique_sequences.txt", 'w') as file_out:
    # Ouvre le fichier d'entrée
    with open(filename,"r") as file:

        # Écrit les lignes d'en-tête dans le fichier de sortie
        for line in file:
            if line.startswith("BLASTN") or line.startswith("Reference") or line.startswith("Database") or line.startswith("Query=") or line.startswith("Length="):
                file_out.write(line)

        # Ajoute un saut de ligne après les lignes d'en-tête
        file_out.write("\n")

        # Parcourt toutes les lignes du fichier
        for i in range(0, len(lines_list)-1):

            # Parcourt toutes les requêtes
            for x in range(1, int(num)+1):

                # Extrait les lignes correspondant à la requête actuelle
                lineA = lines_list[i]
                lineB = lines_list[i+1]

                # Si la ligne A contient la requête actuelle et que la ligne B est un saut de ligne, ajoute la ligne A aux nouvelles données
                if f'Query_{x}' in lineA and '\n' == lineB :
                    new_data.append(lineA)

        # Écrit les nouvelles données dans le fichier de sortie
        for line in new_data:
            file_out.write(line)

# Ouvre le fichier de sortie pour ajouter des sauts de ligne
with open("unique_sequences.txt", 'r') as file:
    lines_list = file.readlines()

new_data = []
start_adding_newlines = False

# Parcourt toutes les lignes du fichier
for i in range(0, len(lines_list)-1):

    # Extrait les lignes actuelles
    lineA = lines_list[i]
    lineB = lines_list[i+1]

    # Commence à ajouter des sauts de ligne à partir de la première occurrence de "Query_1"
    if "Query_1" in lineA:
        start_adding_newlines = True

    # Si on doit ajouter des sauts de ligne
    if start_adding_newlines:
        # Vérifie s'il y a assez d'éléments dans la ligne avant d'essayer d'accéder à l'index 3
        if len(lineA.split()) > 3 and len(lineB.split()) > 1:
            # Extrait le numéro de fin de la ligne A et le numéro de début de la ligne B
            end_num_A = int(lineA.split()[3])
            start_num_B = int(lineB.split()[1])

            # Si le numéro de fin de la ligne A ne suit pas le numéro de début de la ligne B, ajoute un saut de ligne
            if end_num_A + 1 != start_num_B:
                new_data.append(lineA)
                new_data.append('\n')
            else:
                new_data.append(lineA)

# Écrit les nouvelles données dans le fichier de sortie
with open("unique_sequences.txt", 'w') as file_out:
    for line in new_data:
        file_out.write(line)

# Affiche un message indiquant que le script a terminé
print("Done. The unique sequences have been extracted to unique_sequences.txt.")
