#!/bin/bash

# Author: Maxime Naour
# Date(dd-mm-yyyy): 11-10-2023

### Disclaimer : Activer l'environnement contenant contenant les packages edirect, fastx, BLAST+ NCBI et pv avant d'exécuter ce script 
### --> Pour moi, tous les packages nécessaires sont dans mon environnement conda edirect : conda activate edirect
### Ne pas avoir de dossier temp/ dans le dossier courant

# Vérifier si l'environnement conda edirect est activé
if [[ $(conda info --envs | grep -c "^\s*edirect\s") -eq 0 ]]; then
  echo "L'environnement conda edirect (composé des packages edirecte, fastx et BLAST+ NCBI) doit être activé pour exécuter ce script."
  echo -e "Veuillez activer l'environnement conda edirect pour utiliser les commandes edirect, fastx et ncbi (BLAST) en utilisant la commande : $(tput setaf 1)conda activate edirect$(tput sgr0)\n"
  echo "~~> edirect : $(tput setaf 4)conda install -n edirect -c bioconda entrez-direct$(tput sgr0)"
  echo "~~> fastx toolkit : $(tput setaf 4)conda install -n edirect -c bioconda fastx_toolkit$(tput sgr0)"
  echo "~~> ncbi(BLAST Command Line) dans l'environnement edirect : $(tput setaf 4)conda install -n edirect -c bioconda blast$(tput sgr0)"
  echo "~~> pv dans l'environnement edirect : $(tput setaf 4)conda install -n edirect -c conda-forge pv$(tput sgr0)"
  exit 1
fi

#**********************************************************
# Partie à modifier en fonction de l'espèce d'intérêt
#**********************************************************

# Tableau contenant les espèces d'intérêt
# Exemple : species=("Rattus norvegicus" "Homo sapiens" "Mus musculus")

species=("Rattus norvegicus" "Homo sapiens" "Mus musculus" "viral")

# Si vous voulez télécharger toutes les séquences RefSeq virales identifiées à ce jour, ajoutez dans le tableau "species" ci-dessus "viral".

#**********************************************************
# END * END * END * END * END * END * END * END * END * END
#**********************************************************

### Début du script ###

echo "*************************"
echo "**** blast-genome.sh ****"
echo "*created by Maxime Naour*"
echo "****** version 1.0 ******"
echo "******  12/10/2023 ******"
echo "*************************"
echo ""
echo "Pour toutes questions, vous pouvez me contacter aux adresses mails suivantes : maxime.naour@inrae.fr *#* naour.maxime@gmail.com"
echo ""
echo "Ce script automatisé permet dans l'ordre de : "
echo "1) Télécharger des génomes RefSeq (indiquez les espèces d'intérêt à la ligne 28 du script)"
echo "2) Créer une base de données à partir de génomes déjà présents ou nouvellement téléchargés"
echo "3) Transformer un fichier FASTQ d'intérêt en FASTA pour réaliser un BLAST sur le génome des espèces choisies"
echo ""

# Créer le dossier temp/ s'il n'existe pas déjà
if [ ! -d "temp/" ]; then
  mkdir temp/
fi

# Initialiser la variable pour savoir si au moins un génome a été téléchargé
downloaded=false

# Demander à l'utilisateur s'il a besoin de télécharger un ou plusieurs génomes
read -p "Avez-vous besoin de télécharger un ou plusieurs génomes ? (y/n) " download_genomes

case $download_genomes in
  [Yy]* )
    # Télécharger les génomes RefSeq de toutes les espèces d'intérêt
    for specie in "${species[@]}"
    do
      # Afficher le nom de l'espèce d'intérêt
      echo -e "Espèce d'intérêt : \033[31m$specie\033[0m"

      # Trouver les liens https du génome RefSeq de l'espèce d'intérêt
      if [ "$specie" = "viral" ]; then
        links="https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz"
      else
        links=$(esearch -db assembly -query "$specie[ORGN]" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F'/' '{print $0"/"$NF"_genomic.fna.gz"}' | sed 's|ftp://|https://|')
      fi

      # Afficher les génomes proposés par la NCBI pour l'espèce d'intérêt
      echo -e "Voici les génomes proposés par la NCBI pour \033[31m$specie\033[0m :"
      echo "$links"

      # Télécharger le génome si nécessaire
      while true; do
        if [ "$specie" = "viral" ]; then
          read -p $'Avez-vous besoin de télécharger le génome viral ? (y/n) ' download
        else
          read -p $'Avez-vous besoin de télécharger un des génomes de \e[31m'"$specie"$'\e[0m ? (y/n) ' download
        fi

        case $download in
          [Yy]* )
            # Laisser l'utilisateur choisir le génome à télécharger
            if [ "$specie" = "viral" ]; then
              genome_link="https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz"
            else
              read -p "Entrez le lien https du génome que vous souhaitez télécharger : " genome_link
            fi

            # Extraire le nom du fichier à partir du lien
            genome_file=$(basename "$genome_link")

            # Télécharger le génome choisi par l'utilisateur dans le dossier temp/ en utilisant axel pour le téléchargement en parallèle
            echo "Téléchargement du génome en cours..."
            axel -n 10 -o temp/ "$genome_link"

            # Mettre la variable à true pour indiquer qu'au moins un génome a été téléchargé
            downloaded=true

            break
            ;;
          [Nn]* )
            # Passer à l'espèce suivante si l'utilisateur ne veut pas télécharger le génome
            if [ "$downloaded" = true ]; then
              echo "Passage à l'espèce suivante..."
            fi
            break
            ;;
          * )
            if [ "$specie" = "viral" ]; then
              echo "Veuillez répondre par 'y' ou 'n'."
            else
              echo "Veuillez répondre par 'y' ou 'n'. (Ctrl + c pour quitter)"
            fi
            ;;
        esac
      done
    done

    # Demander à l'utilisateur s'il souhaite concaténer les fichiers génomiques
    while true; do
      read -p "Voulez-vous concaténer les fichiers génomiques ? (y/n) " concat_genomes
      case $concat_genomes in
        [Yy]* )
          # Concaténer tous les fichiers dans un seul fichier si au moins un génome a été téléchargé
          if [ "$downloaded" = true ]; then
            # Demander à l'utilisateur combien de fichiers il souhaite concaténer
            read -p "Combien de fichiers souhaitez-vous concaténer pour constituer la base de données ? " num_files

            for ((i=1; i<=$num_files; i++))
            do
              read -p "Nom du fichier $i (avec extension \".fna.gz\") : " file
              # Vérifier si le fichier existe dans le dossier temp/
              if [ ! -f "temp/$file" ]; then
                echo "Le fichier $file n'existe pas dans le dossier temp/."
                exit 1
              fi
              files="$files temp/$file"
            done

            # Concaténer les fichiers sélectionnés dans un seul fichier
            echo "Concaténation des fichiers en cours..."
            cat $files > temp/genomes_DB.fna.gz

            # Décompresser le fichier contenant tous les génomes
            echo "Décompression du fichier contenant tous les génomes en cours..."
            gunzip -k temp/genomes_DB.fna.gz
          else
            echo "Aucun génome n'a été téléchargé. Vous devez fournir un fichier d'entrée nommé \"genomes_DB.fna\" (fichier décompressé) dans le dossier temp/ pour créer la base de données génomique."
          fi
          break
          ;;
        [Nn]* )
          # Vérifier si le fichier "genomes_DB.fna" se trouve dans le répertoire temp/
          if [ -f "temp/genomes_DB.fna" ]; then
            break
          else
            echo -e "\033[31mATTENTION : Vous devez fournir un fichier d'entrée nommé \"genomes_DB.fna\" (fichier décompressé) dans le dossier temp/ pour créer la base de données génomique. \nFin du processus.\033[0m"
            exit 1
          fi
          ;;
        * )
          echo "Veuillez répondre par 'y' ou 'n'."
          ;;
      esac
    done
    ;;
  [Nn]* )
    echo "Si vous avez déjà téléchargé le ou les génomes d'intérêt, assurez-vous qu'il se trouve dans le dossier temp/ nouvellement créé (répertoire où se trouve le script)."
    # Afficher les fichiers présents dans le dossier temp/ 
    echo -e "\nVoici les fichiers présents dans le dossier temp/ :"
    ls temp/*.fna.gz | xargs -n 1 basename | tr '\n' ' ' | sed 's/ /   /g' && echo

    # même affichage avec awk
    # ls temp/*.fna.gz | awk -F'/' '{print $NF}' | tr '\n' ' ' | sed 's/ /   /g' && echo
    
    echo ""
    # Demander à l'utilisateur s'il souhaite concaténer les fichiers génomiques
    while true; do
      read -p "Voulez-vous concaténer les fichiers génomiques ? (y/n) " concat_genomes
      case $concat_genomes in
        [Yy]* )
          # Demander à l'utilisateur combien de fichiers il souhaite concaténer
          read -p "Combien de fichiers souhaitez-vous concaténer pour constituer la base de données ? " num_files

          for ((i=1; i<=$num_files; i++))
          do
            read -p "Nom du fichier $i (avec extension \".fna.gz\") : " file
            # Vérifier si le fichier existe dans le dossier temp/
            if [ ! -f "temp/$file" ]; then
              echo "Le fichier $file n'existe pas dans le dossier temp/."
              exit 1
            fi
            files="$files temp/$file"
          done

          # Concaténer les fichiers sélectionnés dans un seul fichier
          echo "Concaténation des fichiers en cours..."
          cat $files > temp/genomes_DB.fna.gz

          # Décompresser le fichier contenant tous les génomes
          echo "Décompression du fichier contenant tous les génomes en cours..."
          pv temp/genomes_DB.fna.gz | gunzip -k > temp/genomes_DB.fna
          break
          ;; 
        [Nn]* )
          # Vérifier si le fichier "genomes_DB.fna" se trouve dans le répertoire temp/
          if [ -f "temp/genomes_DB.fna" ]; then
            break
          else
            echo -e "\033[31mATTENTION : Vous devez fournir un fichier d'entrée nommé \"genomes_DB.fna\" (fichier décompressé) dans le dossier temp/ pour créer la base de données génomique. \nFin du processus.\033[0m"
            exit 1
          fi
          ;;
        * )
          echo "Veuillez répondre par 'y' ou 'n'."
          ;;
      esac
    done
    ;;
  * )
    echo "Veuillez répondre par 'y' ou 'n'."
    ;;
esac

# Créer un fichier d'index pour la base de données de BLAST
while true; do
  read -p "Voulez-vous créer un fichier d'index pour la base de données de BLAST ? (y/n) " create_index
  case $create_index in
    [Yy]* )
      echo "Création d'un fichier d'index pour la base de données de BLAST en cours..."
      makeblastdb -in temp/genomes_DB.fna -dbtype nucl -out temp/genomes_DB
      break
      ;;
    [Nn]* )
      break
      ;;
    * )
      echo "Veuillez répondre par 'y' ou 'n'. (Ctrl + c pour quitter)"
      ;;
  esac
done

# Demander à l'utilisateur le nom du fichier FASTQ pour le query du BLAST
while true; do
  read -p $'Entrez le nom du fichier FASTQ pour le \033[0;31mquery\033[0m du BLAST (avec l\'extension .fastq ou .fq) : ' fastq_file

  # Vérifier si le fichier FASTQ existe
  if [ ! -f "$fastq_file" ]; then
    echo "Le fichier $fastq_file n'existe pas. Veuillez entrer un nom de fichier valide. (Ctrl + c pour quitter)"
  else
    break
  fi
done

# Vérifier le format du fichier FASTQ
while true; do
  read -p $'Voulez-vous vérifier le format du fichier \033[0;31mFASTQ\033[0m ? (y/n) ' check_fastq
  case $check_fastq in
    [Yy]* )
      echo -e "Vérification du format du fichier \033[0;31mFASTQ\033[0m en cours..."
      fastq_quality_filter -Q33 -v -i $fastq_file -o /dev/null
      break
      ;;
    [Nn]* )
      break
      ;;
    * )
      echo "Veuillez répondre par 'y' ou 'n'. (Ctrl + c pour quitter)"
      ;;
  esac
done

# Convertir le fichier FASTQ en FASTA
while true; do
  read -p $'Voulez-vous convertir le fichier \033[0;31mFASTQ\033[0m en \033[0;31mFASTA\033[0m ? (y/n) ' convert_fastq
  case $convert_fastq in
    [Yy]* )
      echo -e "Conversion du fichier \033[31mFASTQ\033[0m en \033[31mFASTA\033[0m en cours..."
      awk 'NR%4==1 {print ">"$1} NR%4==2 {print $1}' $fastq_file > ${fastq_file%.*}.fna
      break
      ;;
    [Nn]* )
      break
      ;;
    * )
      echo "Veuillez répondre par 'y' ou 'n'. (Ctrl + c pour quitter)"
      ;;
  esac
done

# Réaliser le BLAST sur le fichier FASTA avec la base de données précédemment créée
while true; do
  read -p $'Voulez-vous réaliser le BLAST ? (y/n) ' run_blast
  case $run_blast in
    [Yy]* )
      # Check if the BLAST database is valid
      if ! blastdbcmd -info -db temp/genomes_DB >/dev/null 2>&1; then
        echo "La base de données BLAST n'est pas valide."
        exit 1
      fi

      # Determine the number of batches needed
      total_sequences=$(grep -c "^>" "${fastq_file%.*}.fna")
      num_batches=$((($total_sequences + 9999) / 10000))

      # Create a temporary directory for batch files
      temp_dir=$(mktemp -d)
      echo "Dossier temporaire créé : $temp_dir"
      echo "Un fichier log est disponible dans le répertoire courant pour voir les erreurs potentiels lors du BLAST : blast_log.txt"

      echo "Nombre total de fichiers batch : $num_batches"

      # Split the input file into batches
      for ((i=0; i<$num_batches; i++)); do
          start=$((i * 10000 + 1))
          end=$((start + 9999 - 1))
          batch_file="$temp_dir/batch_$i.fasta"
          sed -n "$start,$end p" "${fastq_file%.*}.fna" > "$batch_file"

          # Check if the last line of the batch file starts with ">@"
          last_line=$(tail -n 1 "$batch_file")
          if [[ $last_line == ">@"* ]]; then
              # Move the last line to the next batch file
              next_batch_file="$temp_dir/batch_$((i+1)).fasta"
              echo "$last_line" >> "$next_batch_file"
              sed -i '$d' "$batch_file"
          fi

          echo "Batch $((i+1))/$num_batches en cours de création."
      done

      # Run BLAST on each batch file in parallel
      blast_output="blast_output.txt"
      batch_count=0
      for batch_file in "$temp_dir"/batch_*.fasta; do
          ((batch_count++))
          echo "Batch $batch_count/$num_batches en cours de BLAST. Cette procédure peut durer plusieurs minutes..."
          sequence_count=1
          while IFS= read -r line; do
              if [[ $line == ">"* ]]; then
                  ((sequence_count++))
              fi
          done < "$batch_file"

          echo -e "Lancement de BLAST pour le Batch $batch_count/$num_batches en utilisant $(nproc) coeurs."
          blastn -query "$batch_file" -db temp/genomes_DB -out "$temp_dir/blast_output_${batch_count}.txt" -num_threads "$(nproc)" -num_descriptions 5 -num_alignments 5 2>&1 | tee -a blast_log.txt
          if [ "${PIPESTATUS[0]}" -ne 0 ]; then
              echo "Erreur lors de l'exécution de BLAST pour le Batch $batch_count/$num_batches"
              cat "$temp_dir/blast_output_${batch_count}_log.txt"
          else
              echo "BLAST du Batch $batch_count/$num_batches terminé."
          fi
      done

      # Concatenate the BLAST output files
      cat "$temp_dir"/blast_output_*.txt > "$blast_output"

      # Process the concatenated output file
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
      done < "$blast_output" > output_space.txt

      # Add the BLAST summary to the output file
      echo "~ Summary of BLAST ~" > temp_output.txt
      echo "Total number of blasted sequences : $((sequence_count-1))" >> temp_output.txt
      echo "Number of sequences with a hit : $hit_count" >> temp_output.txt
      echo "Number of sequences without hits : $no_hit_count" >> temp_output.txt
      echo "~ ~ ~ ~ ~ ~ ~ ~ ~ ~" >> temp_output.txt
      echo "" >> temp_output.txt
      if [ -s output_space.txt ]; then
        cat output_space.txt >> temp_output.txt

        # Définir le nom du fichier
        input_file="temp_output.txt"

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

        # Split the output file into multiple files of 1 GB
        if [ "$(stat -c%s temp_output.txt)" -gt 1000000000 ]; then
          split -b 1G --additional-suffix=.txt --numeric-suffixes temp_output.txt "processed_${fastq_file%.*}_part"
          num_parts=$(ls -1 processed_"${fastq_file%.*}"_part* | wc -l)
          for ((i=1; i<=$num_parts; i++)); do
              part_file="processed_${fastq_file%.*}_part${i}.txt"
              mv "$part_file" "processed_${fastq_file%.*}_part${i}.txt"
          done
          echo "Le fichier de sortie a été divisé en plusieurs fichiers de 1 Go maximum."
          echo "Les fichiers de sortie ont été renommés en processed_${fastq_file%.*}_partX.txt et sont disponibles dans le répertoire courant."
        else
          mv temp_output.txt "processed_${fastq_file%.*}.txt"
          echo -e "\nLe fichier de sortie a été renommé en \033[0;34mprocessed_${fastq_file%.*}.txt\033[0m et est disponible dans le répertoire courant."
          echo -e "\033[0;34mLe processus est terminé.\033[0m"
          echo -e "\033[0;34mMerci d'avoir utilisé le script blast-genome.sh.\033[0m"
        fi
      else
        echo "Le fichier de sortie est vide."
      fi

      # Clean up intermediate files
      rm -rf "$temp_dir"

      break
      ;;
    [Nn]* )
      echo -e "\033[0;34mLe processus est terminé.\033[0m"
      echo -e "\033[0;34mMerci d'avoir utilisé le script blast-genome.sh.\033[0m"
      break
      ;;
    * )
      echo "Veuillez répondre par 'y' ou 'n'. (Ctrl + c pour quitter)"
      ;;
  esac
done
