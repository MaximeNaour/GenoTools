#!/bin/bash

# Automatisation non fonctionnelle. Il réalise la comparaison uniquement pour le premier couple de numéro d'accès. 

input_file="input.csv"
output_file="compare_outputs.txt"

# Supprime le fichier de sortie s'il existe déjà
if [ -f "$output_file" ]; then
  rm "$output_file"
fi

# Convertir les fins de ligne du fichier d'entrée en format Unix
dos2unix $input_file

# Parcourez chaque ligne du fichier d'entrée
count=0
tail -n +2 "$input_file" | while IFS=',' read -r refseq_assembly ncbi_reference_sequence
do
  count=$((count+1))
  echo "Traitement de la ligne $count : $refseq_assembly, $ncbi_reference_sequence"
  
  # Récupérer l'URL du fichier assembly_summary.txt pour l'organisme d'intérêt
  summary_url=$(curl -s "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt" -o - | awk -F "\t" -v acc="$refseq_assembly" '$1 == acc {print $20}')

  # Récupérer le nom du fichier FASTA à partir de l'URL du résumé
  fasta_file_name=$(basename ${summary_url}_genomic.fna.gz)

  # Télécharger le fichier FASTA
  curl -O "${summary_url}/${fasta_file_name}" </dev/null

  # Décompresser le fichier FASTA
  gzip -d ${fasta_file_name}

  # Lire le fichier FASTA décompressé
  fasta_file_content=$(cat ${fasta_file_name%.gz})

  # Récupérer les séquences pour les deux identifiants de séquence
  sequence1=$(echo "$fasta_file_content")
  sequence2=$(esearch -db nuccore -query "$ncbi_reference_sequence" | efetch -format fasta)

  # Supprimer les informations d'en-tête des séquences
  sequence1_clean=$(echo "$sequence1" | grep -v ">" | tr -d '\n' | tr -d '\r')
  sequence2_clean=$(echo "$sequence2" | grep -v ">" | tr -d '\n' | tr -d '\r')

    # Compter le nombre de nucléotides pour chaque séquence
  sequence1_length=${#sequence1_clean}
  sequence2_length=${#sequence2_clean}

  # Comparer les nombres de nucléotides et écrire le résultat dans le fichier de sortie
  echo "Comparaison pour $refseq_assembly et $ncbi_reference_sequence" >> "$output_file"
  if [ "$sequence1_length" -eq "$sequence2_length" ]; then
    echo "Le nombre total de nucléotides est identique pour $refseq_assembly et $ncbi_reference_sequence." >> "$output_file"
    echo "$refseq_assembly : $sequence1_length nucléotides" >> "$output_file"
    echo "$ncbi_reference_sequence : $sequence2_length nucléotides" >> "$output_file"
  else
    echo "Le nombre total de nucléotides n'est pas identique pour $refseq_assembly et $ncbi_reference_sequence." >> "$output_file"
    echo "$refseq_assembly : $sequence1_length nucléotides" >> "$output_file"
    echo "$ncbi_reference_sequence : $sequence2_length nucléotides" >> "$output_file"
  fi

  # Ajouter un saut de ligne dans le fichier de sortie pour séparer chaque comparaison
  echo "" >> "$output_file"

  # Supprimer le fichier FASTA téléchargé
  rm ${fasta_file_name%.gz}

done
