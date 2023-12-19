# Protocole d'Identification de Primers Spécifiques

Ce protocole détaille les étapes pour identifier des primers spécifiques à l'échelle de la souche bactérienne (Protocole 1) et à l'échelle de l'espèce bactérienne (Protocole 2). 

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

5. **Trouve tous les fichiers .fna.gz qui contiennent le nom de la souche d'intérêt (plusieurs motifs) dans leur première ligne, affiche leur nom et la première ligne de chaque fichier, et compte le nombre de ces fichiers**
```
{ find gbct_[nom_classe]/ -name "*.fna.gz" -print0 | xargs -0 -I{} bash -c 'if gunzip -c "{}" | head -1 | grep -q "[motif1]\|[motif2]\|[motif3]"; then echo "{}"; gunzip -c "{}" | head -1; fi'; } | tee >(grep -v '^>' | wc -l | xargs -I{} echo "Nombre de fichiers = {}")
```

6. **Supprimer le ou les fichiers FASTA identifiés appartenant à la souche d'intérêt**
```
rm -rf gbct_[nom_classe]/[fichier1] gbct_[nom_classe]/[fichier2] (etc...)
```

### Alternative 

#### Utilisation du catalogue Mouse Gastrointestinal Bacteria Catalogue (MGBC) si le génome de cette souche provient de ce catalogue

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
mkdir -p [nom_espèce]_primers/gbct_[nom_espèce]
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

4. **Répertoire pour les génomes qui seront inclus dans la base de données**
```
mkdir gbct_[nom_classe]
```

### Etape 2 : Téléchargement des génomes
0. **Activer l'environnement conda contenant l'outil Entrez-Direct** (path = /usr/local/genome/Anaconda3/envs/entrez-direct-15.6)
```
conda activate entrez-direct-15.6
```

1. **Rechercher la souche d'intérêt dans la base de données RefSeq de NCBI et récupérer les liens FTP des génomes RefSeq**
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

5. **Trouve tous les fichiers .fna.gz qui contiennent le nom de la souche d'intérêt (plusieurs motifs) dans leur première ligne, affiche leur nom et la première ligne de chaque fichier, et compte le nombre de ces fichiers**
```
{ find gbct_[nom_classe]/ -name "*.fna.gz" -print0 | xargs -0 -I{} bash -c 'if gunzip -c "{}" | head -1 | grep -q "[motif1]\|[motif2]\|[motif3]"; then echo "{}"; gunzip -c "{}" | head -1; fi'; } | tee >(grep -v '^>' | wc -l | xargs -I{} echo "Nombre de fichiers = {}")
```

6. **Supprimer le ou les fichiers FASTA identifiés appartenant à la souche d'intérêt**
```
rm -rf gbct_[nom_classe]/[fichier1] gbct_[nom_classe]/[fichier2] (etc...)
```

### Alternative 

#### Utilisation du catalogue Mouse Gastrointestinal Bacteria Catalogue (MGBC) si le génome de cette souche provient de ce catalogue

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

4. **Activer l'environnement conda ayant l'outil makeblastdb** (ex : path = /usr/local/genome/Anaconda3/envs/blast-2.13.0)
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

### Étape finale : Conception des amorces/ primers
- Examiner les fichiers de sortie pour identifier des primers potentiels à l'aide de l'outil "PCR Primer Design" (Eurofins) : https://eurofinsgenomics.eu/en/ecom/tools/pcr-primer-design/

Voici les paramètres et valeurs fréquemment utilisés (à adapter selon votre situation) : 

![image](https://github.com/MaximeNaour/GenoTools/blob/main/images/criteria.png?raw=true)
