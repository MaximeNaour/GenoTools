#!/bin/bash

# Author: Maxime Naour
# Date(dd-mm-yyyy): 11-10-2023

# Disclaimer : Activer les environnements conda edirect et ncbi : conda activate edirect ncbi

# Vérifier si l'environnement conda edirect est activé
if [[ $(conda info --envs | grep -c "^\s*edirect\s") -eq 0 ]]; then
  echo "L'environnement conda edirect (composé des packages edirecte, fastx et BLAST+ NCBI) doit être activé pour exécuter ce script."
  echo -e "Veuillez activer l'environnement conda edirect pour utiliser les commandes edirect, fastx et ncbi (BLAST) en utilisant la commande : $(tput setaf 1)conda activate edirect$(tput sgr0)\n"
  echo "~~> edirect : $(tput setaf 4)conda install -n edirect -c bioconda entrez-direct$(tput sgr0)"
  echo "~~> fastx toolkit : $(tput setaf 4)conda install -n edirect -c bioconda fastx_toolkit$(tput sgr0)"
  echo "~~> ncbi(BLAST Command Line) dans l'environnement edirect : $(tput setaf 4)conda install -n edirect -c bioconda blast$(tput sgr0)"
  exit 1
fi

#*********************************************************
# Partie à modifier en fonction de l'espèce d'intérêt
#*********************************************************

# Nom de l'espèce d'intérêt
specie="Rattus norvegicus"

#*********************************************************
# END * END * END * END * END * END * END * END * END * END
#*********************************************************

### Début du script ###

# Afficher le nom de l'espèce d'intérêt de manière plus lisible
echo "Espèce d'intérêt : $specie"

# Trouver les liens https du génome complet RefSeq de l'espèce d'intérêt
links=$(esearch -db assembly -query "$specie[ORGN]" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F'/' '{print $0"/"$NF"_genomic.fna.gz"}' | sed 's|ftp://|https://|')

# Afficher les génomes complets proposés par la NCBI pour l'espèce d'intérêt
echo "Voici les génomes complets proposés par la NCBI pour $specie :"
echo "$links"

# Demander à l'utilisateur s'il a besoin de télécharger le génome complet
read -p "Avez-vous besoin de télécharger le génome complet ? (y/n) " download

if [ "$download" = "y" ]; then
  # Laisser l'utilisateur choisir le génome complet à télécharger
  read -p "Entrez le lien https du génome complet que vous souhaitez télécharger : " genome_link

  # Télécharger le génome complet choisi par l'utilisateur dans le dossier temp/ en utilisant axel pour le téléchargement en parallèle
  echo "Téléchargement du génome complet en cours..."
  axel -n 10 -o temp/ $genome_link
fi

# Lister les fichiers du dossier temp/ pour que l'utilisateur choisisse le fichier de référence pour le BLAST
echo "Voici les fichiers présents dans le dossier temp/ :"
ls temp/

# Demander à l'utilisateur le nom du fichier de référence pour le BLAST
read -p "Entrez le nom du fichier de référence pour le BLAST (avec l'extension .fna.gz) : " filename

# Demander à l'utilisateur le nom du fichier FASTQ pour le query du BLAST
read -p "Entrez le nom du fichier FASTQ pour le query du BLAST (avec l'extension .fastq ou .fq) : " fastq_file 

# Vérifier le format du fichier FASTQ
read -p "Voulez-vous vérifier le format du fichier FASTQ ? (y/n) " check_fastq
if [ "$check_fastq" = "y" ]; then
  echo "Vérification du format du fichier FASTQ en cours..."
  fastq_quality_filter -Q33 -v -i $fastq_file -o /dev/null
fi

# Convertir le fichier FASTQ en FASTA
read -p "Voulez-vous convertir le fichier FASTQ en FASTA ? (y/n) " convert_fastq
if [ "$convert_fastq" = "y" ]; then
  echo "Conversion du fichier FASTQ en FASTA en cours..."
  awk 'NR%4==1 {print ">"$1} NR%4==2 {print $1}' $fastq_file > ${fastq_file%.*}.fna
fi

# Décompresser le fichier de référence pour le BLAST
read -p "Voulez-vous décompresser le fichier de référence pour le BLAST ? (y/n) " decompress_ref
if [ "$decompress_ref" = "y" ]; then
  echo "Décompression du fichier de référence pour le BLAST en cours..."
  gunzip -k temp/$filename
fi

# Créer un fichier d'index pour la base de données de BLAST
read -p "Voulez-vous créer un fichier d'index pour la base de données de BLAST ? (y/n) " create_index
if [ "$create_index" = "y" ]; then
  echo "Création d'un fichier d'index pour la base de données de BLAST en cours..."
  makeblastdb -in temp/${filename%.gz} -dbtype nucl -out temp/${filename%.*.gz}
fi

# Réaliser le BLAST avec les outils BLAST de la NCBI en utilisant le multithreading
read -p "Voulez-vous réaliser le BLAST ? (y/n) " run_blast
if [ "$run_blast" = "y" ]; then
  echo "Réalisation du BLAST en cours..."
  blastn -query ${fastq_file%.*}.fna -db temp/${filename%.*.gz} -out output.txt -num_threads 6 

  # Modifier le fichier de sortie pour afficher les résultats pour chaque séquence
  echo "Modification du fichier de sortie en cours..."
  sequence_count=1
  hit_count=0
  no_hit_count=0
  while IFS= read -r line; do
      if [[ $line == "Query="* ]]; then
          echo "*****************************"
          echo "Result sequence $sequence_count"
          echo "*****************************"
          echo "$line"
          ((sequence_count++))
      elif [[ $line == "***** No hits found *****" ]]; then
          echo "$line"
          ((no_hit_count++))
      elif [[ $line == "Sequences producing significant alignments:"* ]]; then
          echo "$line"
          ((hit_count++))
      else
          echo "$line"
      fi
  done < "output.txt" > output_space.txt

  # Ajouter les informations générales du BLAST en haut du fichier
  echo "~ Summary of BLAST ~" > temp_output.txt
  echo "Total number of blasted sequences : $((sequence_count-1))" >> temp_output.txt
  echo "Number of sequences with a hit : $hit_count" >> temp_output.txt
  echo "Number of sequences without hits : $no_hit_count" >> temp_output.txt
  echo "~ ~ ~ ~ ~ ~ ~ ~ ~ ~" >> temp_output.txt
  echo "" >> temp_output.txt
  cat output_space.txt >> temp_output.txt
  mv temp_output.txt processed_${fastq_file%.*}.txt && echo "Le fichier de sortie a été renommé en processed_${fastq_file%.*}.txt"
  # Supprimer les fichiers output.txt et output_space.txt
  rm output.txt output_space.txt
fi

# Définir le nom du fichier
input_file="processed_${fastq_file%.*}.txt"

# Extraire le bloc de texte à déplacer
block_to_move=$(awk '/~ Summary of BLAST ~/,/~ ~ ~ ~ ~ ~ ~ ~ ~ ~/ {print}' "$input_file")

# Supprimer les six premières lignes du fichier
sed -i '1,6d' "$input_file"

# Définir un motif partiel
partial_pattern="Result sequence 1"

# Trouver la ligne correspondant au motif partiel
line_number=$(grep -n "$partial_pattern" "$input_file" | cut -d: -f1)

# Vérifier si la ligne a été trouvée avant d'effectuer l'insertion
if [ -n "$line_number" ]; then
    # Calculer le numéro de la ligne où insérer le bloc
    insert_line=$((line_number - 2))
    
    # Insérer le bloc de texte deux lignes avant la ligne trouvée
    awk -v block="$block_to_move" -v line="$insert_line" 'NR==line {print block} {print}' "$input_file" > temp_file.txt && mv temp_file.txt "$input_file"
else
    echo "Pattern non trouvé dans le fichier."
fi

echo "Le processus est terminé."

