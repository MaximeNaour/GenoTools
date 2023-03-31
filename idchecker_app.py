#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
Created on 31/03/2023
@author: maxime.naour@inrae.fr
"""

# pip install biopython

# install application : pyinstaller --onefile --noconsole IdChecker.py

import csv
import os
from tkinter import *
from tkinter import filedialog
from tkinter import ttk, messagebox
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

def split_species_strain(ncbi_name):
    split_name = ncbi_name.split()
    species = split_name[0] + ' ' + split_name[1]
    strain = ' '.join(split_name[2:])
    return species, strain

def update_bacteria_names(input_file, output_file, progress_window):
    with open(input_file, "r") as input_csv, open(output_file, "w", newline='') as output_csv:
        reader = csv.reader(input_csv, delimiter=' ')
        writer = csv.writer(output_csv, delimiter=',')
        writer.writerow(['Specie', 'Strain'])

        num_rows = sum(1 for row in reader)
        input_csv.seek(0)

        progress_var = DoubleVar()
        progress_bar = ttk.Progressbar(progress_window, variable=progress_var, length=200, mode='determinate')
        progress_bar.pack(padx=10, pady=10)

        for i, row in enumerate(reader):
            species = row[0] + ' ' + row[1]
            strain_info = ' '.join(row[2:])

            ncbi_name = search_ncbi(species, strain_info)

            if ncbi_name:
                species, strain = split_species_strain(ncbi_name)
                writer.writerow([species, strain])
            else:
                writer.writerow([species, strain_info])

            progress_var.set(i / num_rows * 100)
            progress_window.update()

        progress_window.destroy()

def browse_file():
    file_path = filedialog.askopenfilename()
    input_file_entry.delete(0, END)
    input_file_entry.insert(END, file_path)

def run_program():
    input_file = input_file_entry.get()

    if not os.path.isfile(input_file):
        messagebox.showerror("Error", "Invalid input file path")
        return

    output_file_path = filedialog.asksaveasfilename(defaultextension='.csv')
    if output_file_path:
        progress_window = Toplevel(root)
        progress_window.title("Processing...")
        progress_window.geometry("250x100")
        progress_window.resizable(False, False)

        progress_label = Label(progress_window, text="Processing...")
        progress_label.pack(padx=10, pady=10)

        update_bacteria_names(input_file, output_file_path, progress_window)
        messagebox.showinfo("Success","Processing completed successfully")

# GUI part
root = Tk()
root.title("Bacteria Names Correction - Developped by Maxime Naour (2023)")
root.geometry("400x150")
root.iconbitmap("images/logo_PhAN.ico")
root.resizable(True, True)

# Frame for input file selection
input_frame = LabelFrame(root, text="Input File")
input_frame.pack(fill=BOTH, expand=True, padx=10, pady=10)

input_file_entry = Entry(input_frame)
input_file_entry.pack(side=LEFT, fill=X, expand=True, padx=5, pady=5)

browse_button = Button(input_frame, text="Browser", command=browse_file)
browse_button.pack(side=LEFT, padx=5, pady=5)

# Run button
run_button = Button(root, text="Run", command=run_program)
run_button.pack(side=BOTTOM, padx=10, pady=10)

# Make widgets resizeable
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)
input_frame.columnconfigure(0, weight=1)
input_file_entry.columnconfigure(0, weight=1)

# Start the application
root.mainloop()
