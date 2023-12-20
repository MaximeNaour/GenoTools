# Protocole d'Identification de Primers Spécifiques

**Réalisé par Maxime Naour**  
**PhD Student - UMR 1280 PhAN (Nantes, France)**  
**date : 20/12/2023**  

**Description : Ce protocole détaille les étapes pour identifier des primers spécifiques à l'échelle de la souche bactérienne (Protocole 1) et à l'échelle de l'espèce bactérienne (Protocole 2).**  

---

## Protocole 1 : Identification de primers spécifiques d'une souche bactérienne

Ce protocole détaille les étapes pour identifier des primers spécifiques à l'échelle de la souche bactérienne.

### Prérequis

- Système ou sous-système Linux (WSL) avec Python 3 installé.
- Environnements Conda avec BLAST+ et Entrez-Direct installés.

### Aide 
#### Identifier tous les environnements conda ayant un outil d'intérêt (par exemple : makeblastdb)
```
conda info --envs | awk '/\/[a-zA-Z0-9]/ {print $2}' | xargs -I {} find {}/bin -name makeblastdb -print
```

### Étape 1: Préparation de l'Environnement de Travail

1. **Créer un Répertoire de Travail et se placer dans ce répertoire**
```
mkdir -p [nom_souche]_primers/gbct_[nom_souche]
cd [nom_souche]_primers
```

2. **Répertoire pour les fichiers de log**
```
mkdir logs
```

3. **Répertoire pour la base de données BLAST**
```
mkdir db
```

4. **Répertoire pour les génomes**
```
mkdir gbct_[nom_classe]
```

### Etape 2 : Téléchargement des génomes
0. **Activer l'environnement conda contenant l'outil Entrez-Direct** (path = /usr/local/genome/Anaconda3/envs/entrez-direct-15.6)
```
conda activate entrez-direct-15.6
```

2. **Rechercher la souche d'intérêt dans la base de données RefSeq de NCBI et récupérer les liens FTP des génomes RefSeq**
```
esearch -db assembly -query "[nom_souche][Organism] AND latest_refseq[filter]" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F'/' '{print $0"/"$NF"_genomic.fna.gz"}' | sed 's|ftp://|https://|'
```

2. **Télécharger le génome de référence de la souche d'intérêt**
```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/[chemin_souche]_genomic.fna.gz -O gbct_[nom_souche]/[nom_souche]_genomic.fna.gz
```

3. **Compter tous les génomes RefSeq de la classe bactérienne de la souche (exemple : Bacteroidia, Clostridia, etc)**
```
esearch -db assembly -query "[nom_classe][Organism] AND latest_refseq[filter]" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F'/' '{print $0"/"$NF"_genomic.fna.gz"}' | sed 's|ftp://|https://|' | wc -l
```

4. **Télécharger tous les génomes RefSeq de la classe bactérienne de la souche (exemple : Bacteroidia, Clostridia, etc)**
```
esearch -db assembly -query "[nom_classe][Organism] AND latest_refseq[filter]" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F'/' '{print $0"/"$NF"_genomic.fna.gz"}' | sed 's|ftp://|https://|' | xargs -n 1 wget -P gbct_[nom_classe]/ 
```
**Ou envoyer cette commande sur un cluster de calcul SGE**
```
qsub -cwd -V -N DL_[nom_classe] -o qlogs.DL_[nom_classe] -e qlogs.DL_[nom_classe] -pe thread 20 -b y "esearch -db assembly -query '[nom_classe][Organism] AND latest_refseq[filter]' | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F'/' '{print $0"/"$NF"genomic.fna.gz"}' | sed 's|ftp://|https://|' | xargs -n 1 wget -P gbct_[nom_classe]/"
```

5. **Trouve tous les fichiers .fna.gz qui contiennent le nom de la souche d'intérêt (plusieurs motifs) dans leur première ligne, affiche leur nom et la première ligne de chaque fichier, et compte le nombre de ces fichiers**
```
{ find gbct_[nom_classe]/ -name "*.fna.gz" -print0 | xargs -0 -I{} bash -c 'if gunzip -c "{}" | head -1 | grep -q "[motif1]\|[motif2]\|[motif3]"; then echo "{}"; gunzip -c "{}" | head -1; fi'; } | tee >(grep -v '^>' | wc -l | xargs -I{} echo "Nombre de fichiers = {}")
```

6. **Supprimer le ou les fichiers FASTA identifiés appartenant à la souche d'intérêt**
```
rm -rf gbct_[nom_classe]/[fichier1].fna.gz gbct_[nom_classe]/[fichier2].fna.gz (etc...)
```

### Alternative 
**Utilisation du catalogue Mouse Gastrointestinal Bacteria Catalogue (MGBC) si le génome de cette souche provient de ce catalogue**

0. **Télécharger le fichier tsv cataloguant tous les MAGS du catalogue MGBC (https://www.sciencedirect.com/science/article/pii/S1931312821005680#mmc4)**  
**Nommer le dit fichier : MGBC_mags.tsv**  
**Le placer dans le répertoire de travail**

1. **Télécharger les génomes des souches spécifiques du catalogue MGBC**
```
qsub -cwd -V -N Download_Genomes -o qlogs -e qlogs -pe thread 30 -b y "mkdir -p [species]_MGBC/ && awk -F '\t' '\$12 == \"Scaffold\" {split(\$1, a, \"_\"); gsub(/\\./, \"\", a[2]); print \"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/\"substr(a[2],1,3)\"/\"substr(a[2],4,3)\"/\"substr(a[2],7,3)\"/\"\$1\"_\"\$2\"/\"\$1\"_\"\$2\"_genomic.fna.gz\"}' MGBC_mags.tsv | xargs -I {} wget -P genomes {}"
```

2. **Compter le nombre de génomes de souches spécifiques téléchargés**
```
ls [species]_MGBC/ | wc -l
```

3. **Trouve tous les fichiers .fna.gz qui contiennent le nom de la souche d'intérêt (plusieurs motifs) dans leur première ligne, affiche leur nom et la première ligne de chaque fichier, et compte le nombre de ces fichiers**
```
{ find genomes/ -name "*.fna.gz" -print0 | xargs -0 -I{} bash -c 'if gunzip -c "{}" | head -1 | grep -q "[motif1]\|[motif2]\|[motif3]"; then echo "{}"; gunzip -c "{}" | head -1; fi'; } | tee >(grep -v '^>' | wc -l | xargs -I{} echo "Nombre de fichiers = {}")
```

4. **Supprimer le ou les fichiers FASTA identifiés appartenant à la souche d'intérêt**
```
rm -rf genomes/[fichier1] (etc...)
```
#### *** Fin de l'alternative ***

### Étape 3 : Préparation de la Base de Données BLAST

1. **Concaténer les génomes téléchargés pour constituer une base de données BLAST**
```
for file in gbct_[nom_classe]/*.fna.gz; do gzip -dc "$file"; done > db/gbct_combined.fna
```
**ou avec les génomes téléchargés à partir du catalogue MGBC**
```
for file in [species]_MGBC/*.fna.gz; do gzip -dc "$file"; done > db/gbct_combined.fna
```
**Ou, Concaténer les génomes téléchargés pour constituer une base de données BLAST sur un cluster SGE**
```
qsub -cwd -V -N Concat_[nom_classe] -o qlogs.Concat_[nom_classe] -e qlogs.Concat_[nom_classe] -pe thread 20 -b y "for file in gbct_[nom_classe]/*.fna.gz; do gzip -dc \$file; done > db/gbct_combined.fna"
```
**ou avec les génomes téléchargés à partir du catalogue MGBC**
```
qsub -cwd -V -N Concat_[species] -o qlogs.Concat_[species] -e qlogs.Concat_[species] -pe thread 20 -b y "for file in [species]_MGBC/*.fna.gz; do gzip -dc \$file; done > db/gbct_combined.fna"
```

2. **Compter le nombre de contigs**
```
rg -c ">" db/gbct_combined.fna
```

3. **Compresser le répertoire contenant les génomes téléchargés pour économiser de l'espace disque puis les supprimer**
```
tar -czvf genomes_[nom_classe].tar.gz gbct_[nom_classe]/ && rm -rf gbct_[nom_classe]/ 
```
**ou avec les génomes téléchargés à partir du catalogue MGBC**
```
tar -czvf [species]_MGBC.tar.gz [species]_MGBC/ && rm -rf [species]_MGBC/
```
**Ou, Compresser le répertoire contenant les génomes téléchargés pour économiser de l'espace disque puis les supprimer sur un cluster SGE**
```
qsub -cwd -V -N Compress_[nom_classe] -o qlogs.Compress_[nom_classe] -e qlogs.Compress_[nom_classe] -pe thread 15 -b y "tar -czvf genomes_[nom_classe].tar.gz gbct_[nom_classe]/ && rm -rf gbct_[nom_classe]/"
```
**ou avec les génomes téléchargés à partir du catalogue MGBC**
```
qsub -cwd -V -N Compress_[species] -o qlogs.Compress_[species] -e qlogs.Compress_[species] -pe thread 15 -b y "tar -czvf [species]_MGBC.tar.gz [species]_MGBC/ && rm -rf [species]_MGBC/"
```

4. **Activer l'environnement conda contenant l'outil makeblastdb** (ex : path = /usr/local/genome/Anaconda3/envs/blast-2.13.0)
```
conda activate blast-2.13.0 
```

5. **Constituer la base de données à partir du fichier multifasta contenant tous les génomes**
```
makeblastdb -in db/gbct_combined.fna -out db/gbct_combined -parse_seqids -dbtype nucl
```

6. **Importer dans le répertoire db/, le génome FASTA décompressé de la bactérie d'intérêt (= Query) afin de faire le BLAST**
```
gzip -dc gbct_[nom_souche]/[nom_souche]_genomic.fna.gz > db/[nom_souche]_genomic.fna
```

7. **Lister les contigs du fichier FASTA de la souche MGBC112867**
```
rg ">" db/[nom_souche]_genomic.fna --no-line-number | awk -F " " '{print $1}' | sed -r 's/>//g' > db/[nom_souche]_contigs.txt
```

### Étape 4 : Réalisation du BLAST

1. **Activer l'environnement conda ayant l'outil blastn** (ex : path = /usr/local/genome/Anaconda3/envs/blast-2.13.0)
```
conda activate blast-2.13.0 
```

2. **Réaliser le BLAST entre le fichier FASTA de la bactérie d'intérêt (= Query) et la base de données précédemment constituées**
```
blastn -db db/gbct_combined -query db/[nom_souche]_genomic.fna -out db/nomatching_[nom_souche].txt -outfmt 2
```
**Ou, execution sur un cluster SGE**
```
qsub -cwd -V -N Blast.[nom_souche] -o qlogs.Blast.[nom_souche] -e qlogs.Blast.[nom_souche] -pe thread 20 -b y "conda activate blast-2.13.0 && blastn -db db/gbct_combined -query db/[nom_souche]_genomic.fna -out db/nomatching_[nom_souche].txt -outfmt 2 && conda deactivate"
```

### Étape 5 : Extraction et Analyse des Séquences Uniques

1. **Lancer le script Python pour identifier les séquences uniques**
```
python unique_seq.py
```

2. **Si vous n'obtenez pas assez de séquences uniques pour votre utilisation, lancer le script Python suivant pour identifier les séquences peu partagées entre votre souche d'intérêt et la base de données**
```
python few_matches.py
```

3. **Conception des amorces/ primers de PCR (cf. étape finale)**  

---

## Protocole 2 : Identification de primers spécifiques à l'échelle de l'espèce

Ce protocole détaille les étapes pour identifier des primers spécifiques à l'échelle de l'espèce bactérienne.

### Prérequis

- Système ou sous-système Linux (WSL) avec Python 3 installé.
- Environnements Conda avec BLAST+, Entrez-Direct, Prokka et Roary installés.

### Étape 1: Préparation de l'Environnement de Travail

1. **Créer un Répertoire de Travail et se placer dans ce répertoire**
```
mkdir [nom_espèce]_primers/
cd [nom_espèce]_primers
```

2. **Répertoire pour les fichiers de log**
```
mkdir logs
```

3. **Répertoire pour la base de données BLAST**
```
mkdir db
```

### Etape 2 : Téléchargement des génomes
0. **Activer l'environnement conda contenant l'outil Entrez-Direct** (path = /usr/local/genome/Anaconda3/envs/entrez-direct-15.6)
```
conda activate entrez-direct-15.6
```

1. **Rechercher l'espèce d'intérêt dans la base de données RefSeq de NCBI et récupérer les liens FTP des génomes RefSeq**
```
esearch -db assembly -query "[nom_espèce][Organism] AND latest_refseq[filter]" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F'/' '{print $0"/"$NF"_genomic.fna.gz"}' | sed 's|ftp://|https://|'
```

2. **Télécharger tous les génomes RefSeq de l'espèce d'intérêt et les placer dans un répertoire**
```
esearch -db assembly -query "[nom_espèce][Organism] AND latest_refseq[filter]" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F'/' '{print $0"/"$NF"_genomic.fna.gz"}' | sed 's|ftp://|https://|' | xargs -n 1 wget -P [nom_espèce]_genomes/ 
```
**Ou envoyer cette commande sur un cluster de calcul SGE**
```
qsub -cwd -V -N DL_[nom_espèce] -o qlogs.DL_[nom_espèce] -e qlogs.DL_[nom_espèce] -pe thread 20 -b y "conda activate entrez-direct-15.6 && esearch -db assembly -query '[nom_espèce][Organism] AND latest_refseq[filter]' | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F'/' '{print \$0\"/\"\$NF\"_genomic.fna.gz\"}' | sed 's|ftp://|https://|' | xargs -n 1 wget -P [nom_espèce]_genomes/ && conda deactivate"
```

3. **Compter le nombre de génomes téléchargés**
```
ls [nom_espèce]_genomes/*.fna.gz | wc -l
```

4. **Compter tous les génomes RefSeq de la classe bactérienne de la souche d'intérêt (exemple : Bacteroidia, Clostridia, etc) sans prendre en compte les génomes RefSeq de l'espèce d'intérêt**
```
esearch -db assembly -query "[nom_classe][Organism] AND latest_refseq[filter] NOT [nom_espèce][Organism]" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F'/' '{print $0"/"$NF"_genomic.fna.gz"}' | sed 's|ftp://|https://|' | wc -l
```

5. **Télécharger tous les génomes RefSeq de la classe bactérienne de la souche d'intérêt et les placer dans un répertoire**
```
esearch -db assembly -query "[nom_classe][Organism] AND latest_refseq[filter] NOT [nom_espèce][Organism]" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F'/' '{print $0"/"$NF"_genomic.fna.gz"}' | sed 's|ftp://|https://|' | xargs -n 1 wget -P gbct_[nom_classe]/ 
```
**Ou envoyer cette commande sur un cluster de calcul SGE**
```
qsub -cwd -V -N DL_[nom_classe] -o qlogs.DL_[nom_classe] -e qlogs.DL_[nom_classe] -pe thread 20 -b y "conda activate entrez-direct-15.6 && esearch -db assembly -query '[nom_classe][Organism] AND latest_refseq[filter] NOT [nom_espèce][Organism]' | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F'/' '{print \$0\"/\"\$NF\"_genomic.fna.gz\"}' | sed 's|ftp://|https://|' | xargs -n 1 wget -P [nom_classe]_genomes/ && conda deactivate"
```

6. **Vérifier que tous les fichiers FASTA ont bien été téléchargés**
```
ls [nom_classe]_genomes/*.fna.gz | wc -l
```

7. **Trouve tous les fichiers .fna.gz qui contiennent l'espèce d'intérêt (si besoin, autres motifs permettant d'identifier les contigs de la bactérie d'intérêt) dans leur première ligne, affiche leur nom et la première ligne de chaque fichier, et compte le nombre de ces fichiers**
```
{ find [nom_classe]_genomes/ -name "*.fna.gz" -print0 | xargs -0 -I{} bash -c 'if gunzip -c "{}" | head -1 | grep -q "[nom_espèce]\|[motif2]\|[motif3]"; then echo "{}"; gunzip -c "{}" | head -1; fi'; } | tee >(grep -v '^>' | wc -l | xargs -I{} echo "Nombre de fichiers = {}")
```
**Normalement, aucun génome de la bactérie d'intérêt doit avoir été identifié mais si un génome est identifié, réaliser la tâche 8.**

8. **Supprimer le ou les fichiers FASTA identifiés appartenant à l'espèce d'intérêt (si besoin)**
```
rm -rf [nom_classe]_genomes/[fichier1].fna.gz [nom_classe]_genomes/[fichier2].fna.gz (etc...)
```

9. **Décompresser temporairement tous les fichiers FASTA du répertoire de classe et les concaténer dans un seul et même fichier multifasta non compressé**
```
qsub -cwd -V -N Concat.FASTA -o qlogs.Concat.FASTA -e qlogs.Concat.FASTA -pe thread 30 -b y "for file in [nom_classe]_genomes/*.fna.gz; do gunzip -c \$file >> db/[nom_classe]_genomes.fna; done"
```

10. **Compter le nombre de contigs dans le fichier multifasta**
```
rg -c ">" db/[nom_classe]_genomes.fna
```

11. **Compresser le répertoire contenant les génomes de la classe bactérienne pour économiser de l'espace disque et le supprimer**
```
qsub -cwd -V -N CompressDirs -o qlogs.TAR -e qlogs.TAR -pe thread 15 -b y "tar -czvf [nom_classe]_genomes.tar.gz [nom_classe]_genomes/ && rm -rf [nom_classe]_genomes/"
```

### Alternative 
**Utilisation du catalogue Mouse Gastrointestinal Bacteria Catalogue (MGBC) si le génome de cette souche provient de ce catalogue (voir Protocole 1)**

### Étape 3 : Identifier le "core genome" de l'espèce d'intérêt
0. **Activer l'environnement conda contenant l'outil Prokka** (path = /usr/local/genome/Anaconda3/envs/entrez-direct-15.6)
```
conda activate prokka-1.14.6
```

1. **Décompresser tous les fichiers FASTA du répertoire "[nom_espèce]_genomes" sur le cluster de calcul**
```
gzip -d [nom_espèce]_genomes/*.fna.gz
```
**Ou envoyer cette commande sur un cluster de calcul SGE**
```
qsub -cwd -V -N DZIP.FASTA -o qlogs.DZIP.FASTA -e qlogs.DZIP.FASTA -pe thread 20 -b y "for file in [nom_espèce]_genomes/*.fna.gz; do gunzip -c \$file > \${file%.gz}; done"
```

2. **Annoter les génomes de l'espèce d'intérêt avec Prokka afin d'obtenir des fichiers GFF3 pour chaque génome**
```
prokka --outdir [nom_espèce]_annotations/[nom_espèce] --prefix [nom_espèce] --genus [genre] --species [espèce] --kingdom Bacteria --gcode 11 --cpus 20 [nom_espèce].fna
```
**Ou envoyer cette commande sur un cluster de calcul SGE**
```
for fna in [nom_espèce]_genomes/*.fna; do qsub -cwd -V -N Prokka_$(basename $fna .fna) -o qlogs.Prokka_$(basename $fna .fna) -e qlogs.Prokka_$(basename $fna .fna) -pe thread 20 -b y "conda activate prokka-1.14.6 && prokka --outdir [nom_espèce]_annotations/$(basename $fna .fna) --prefix $(basename $fna .fna) --genus [genre] --species [espèce] --kingdom Bacteria --gcode 11 --cpus 20 $fna && conda deactivate"; done
```

3. **Activer l'environnement conda contenant l'outil Roary pour l'analyse du "core genome"**
```
conda activate roary-3.13.0
```

4. **Identifier les gènes communs (core genome) et les gènes spécifiques (accessory genome) des génomes de l'espèce d'intérêt**
```
roary -p 20 -f [nom_espèce]_Roary -e --mafft -n -v -cd 95 [nom_espèce]_annotations/*/*.gff
```
**Ou envoyer cette commande sur un cluster de calcul SGE**
```
qsub -cwd -V -N Roary_[nom_espèce] -o qlogs.Roary_[nom_espèce] -e qlogs.Roary_[nom_espèce] -pe thread 40 -b y "conda activate roary-3.13.0 && roary -p 40 -f [nom_espèce]_Roary -e --mafft -n -v -cd 95 [nom_espèce]_annotations/*/*.gff && conda deactivate"
```

5. **Générer un arbre phylogénétique à partir de l'alignement des gènes du core genome**
```
FastTree -nt -gtr [nom_espèce]_Roary/core_gene_alignment.aln > [nom_espèce]_Roary/core_gene_alignment.newick
```
**Ou envoyer cette commande sur un cluster de calcul SGE**
```
qsub -cwd -V -N FastTree_[nom_espèce] -o qlogs.FastTree_[nom_espèce] -e qlogs.FastTree_[nom_espèce] -pe thread 20 -b y "conda activate roary-3.13.0 && FastTree -nt -gtr [nom_espèce]_Roary/core_gene_alignment.aln > [nom_espèce]_Roary/core_gene_alignment.newick && conda deactivate"
```

6. **Générer les graphiques de Roary pour visualiser les résultats**
```
roary_plots.py --format png [nom_espèce]_Roary/core_gene_alignment.newick [nom_espèce]_Roary/gene_presence_absence.csv
```

7. **Lire le fichier summary pour visualiser les résultats de l'analyse pangénomique**
```
more summary.txt
```

8. **Extraire les séquences partagées par tous les génomes de l'espèce d'intérêt (core genome) et placer le fichier FASTA dans le répertoire db/**
```
cp [nom_espèce]_Roary/pan_genome_reference.fa db/[nom_espèce]_core_genome.fa
```

9. **Lister les contigs du "core genome" l'espèce d'intérêt**
```
rg ">" db/[nom_espèce]_core_genome.fa --no-line-number | awk -F " " '{print $1}' | sed -r 's/>//g' > db/[nom_espèce]_contigs.txt
```

10. **Compresser et archiver les répertoires contenant les génomes, annotations et roary pour économiser de l'espace disque puis les supprimer**
```
tar -czvf [nom_espèce].tar.gz [nom_espèce]_annotations/ [nom_espèce]_genomes/ [nom_espèce]_Roary/ && rm -rf [nom_espèce]_annotations/ [nom_espèce]_genomes/
```

### Étape 4 : Préparation de la Base de Données BLAST et Réalisation du BLAST

0. **Activer l'environnement conda contenant l'outil makeblastdb** (ex : path = /usr/local/genome/Anaconda3/envs/blast-2.13.0)
```
conda activate blast-2.13.0 
```

1. **Constituer la base de données BLAST à partir des génomes de la classe bactérienne d'intérêt**
```
makeblastdb -in db/[nom_classe]_genomes.fna -out db/[nom_classe]_genomes -parse_seqids -dbtype nucl
```

2. **Réaliser le BLAST entre le core genome et la base de données BLAST des génomes de la classe d'intérêt**
```
blastn -db db/[nom_classe]_genomes -query db/[nom_espèce]_core_genome.fa -out db/matching_[nom_espèce]_core_genome.txt -outfmt 2
```
**Ou envoyer cette commande sur un cluster de calcul SGE**
```
qsub -cwd -V -N Blast_[nom_espèce]_cg -o qlogs.Blast_[nom_espèce]_cg -e qlogs.Blast_[nom_espèce]_cg -pe thread 20 -b y "conda activate blast-2.13.0 && blastn -db db/[nom_classe]_genomes -query db/[nom_espèce]_core_genome.fa -out db/matching_[nom_espèce]_core_genome.txt -outfmt 2 && conda deactivate"
```

## Étape 7 : Extraction et Analyse des Séquences Uniques pour l'Échelle de l'Espèce

1. **Lancer le script Python pour identifier les séquences uniques spécifiques à l'espèce**
```
python unique_seq.py
```

2. **Si nécessaire, utiliser un script Python pour identifier les séquences peu partagées entre l'espèce d'intérêt et les autres génomes de la classe**
```
python few_matches.py
```

3. **Conception des amorces/primers de PCR pour les séquences uniques identifiées (cf. étape finale)**

---

### Étape finale : Conception des amorces/ primers
- Examiner les fichiers de sortie pour identifier des primers potentiels à l'aide de l'outil "PCR Primer Design" (Eurofins) : https://eurofinsgenomics.eu/en/ecom/tools/pcr-primer-design/

Voici les paramètres et valeurs fréquemment utilisés (à adapter selon votre situation) : 

![image](https://github.com/MaximeNaour/GenoTools/blob/main/images/criteria.png?raw=true)
