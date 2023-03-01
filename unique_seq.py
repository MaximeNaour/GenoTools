#!/save_projet/fine/mnaour/PhAN/Metabarcoding_PhAN/pkg/bin/python3
# -*- coding: utf-8 -*-

"""
Created on 28/02/2023
@author: maxime.naour@inrae.fr
"""

import os

filename = input("Enter your file's name : ")
num = input("Enter the number of queries : ")

lines_list = []
new_data = []

if os.path.exists(filename):

    with open(filename,"r") as file :

        for line in file:
            lines_list = file.readlines()

else:
    print("The document", filename, "is not detected, make sure you place it in your environment")

with open("unique_sequences.txt", 'w') as file_out:
    with open(filename,"r") as file:

        for line in file:
            if line.startswith("BLASTN") or line.startswith("Reference") or line.startswith("Database") or line.startswith("Query=") or line.startswith("Length="):
                file_out.write(line)

        file_out.write("\n")

        for i in range(0, len(lines_list)-1):

            for x in range(1, int(num)+1):

                lineA = lines_list[i]
                lineB = lines_list[i+1]

                if f'Query_{x}' in lineA and '\n' == lineB :

                    new_data.append(lineA)

        for line in new_data:

            file_out.write(line)
            

print("Done. The unique sequences have been extracted to unique_sequences.txt.")