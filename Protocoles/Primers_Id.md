# Protocole d'Identification de Primers Spécifiques

Ce protocole détaille les étapes pour identifier des primers spécifiques à l'échelle de la souche bactérienne (Protocole 1) et à l'échelle de l'espèce bactérienne (Protocole 2). 

## Protocole 1 : Identification de primers spécifiques d'une souche bactérienne

Ce protocole détaille les étapes pour identifier des primers spécifiques à l'échelle de la souche bactérienne.

## Prérequis

- Système ou sous-système Linux (WSL) avec Python 3 installé.
- Environnements Conda avec BLAST+ et Entrez-Direct installés.

## Étape 1: Préparation de l'Environnement de Travail

1. **Créer un Répertoire de Travail et se placer dans ce répertoire**
```
mkdir -p [nom_souche]primers/gbct_[nom_souche]
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

## Etape 2 : Téléchargement des génomes

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

7. **Supprimer le ou les fichiers FASTA identifiés appartenant à la souche d'intérêt**
```
rm -rf gbct_[nom_classe]/[fichier1] gbct_[nom_classe]/[fichier2] (etc...)
```

## Alternative : Utilisation du catalogue Mouse Gastrointestinal Bacteria Catalogue (MGBC) si le génome de cette souche provient de ce catalogue

0. **Télécharger le fichier tsv catalogant tous les MAGS du catalgoue MGBC (https://www.sciencedirect.com/science/article/pii/S1931312821005680#mmc4)**  
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

Étape 3 : Préparation de la Base de Données BLAST

    Concaténer les génomes téléchargés pour constituer une base de données BLAST

    bash

for file in Bacteroidia/*.fna.gz [species]_MGBC/*.fna.gz; do gzip -dc "$file"; done > db/gbct_combined.fna

Création de la base de données BLAST à partir des génomes concaténés

bash

    makeblastdb -in db/gbct_combined.fna -out db/gbct_combined -parse_seqids -dbtype nucl

Étape 4 : Préparation du Génome de Référence pour la Souche d'Intérêt

    Importer le génome de la souche d'intérêt dans le répertoire db

    css

gzip -dc gbct_[nom_souche]/[nom_souche]_genomic.fna.gz > db/[nom_souche]_genomic.fna

Lister les contigs du fichier FASTA de la souche d'intérêt

arduino

    rg ">" db/[nom_souche]_genomic.fna --no-line-number | awk -F " " '{print $1}' | sed -r 's/>//g' > db/[nom_souche]_contigs.txt

Étape 5 : Réalisation du BLAST

    Effectuer le BLAST entre le génome de la souche d'intérêt et la base de données

    bash

blastn -db db/gbct_combined -query db/[nom_souche]_genomic.fna -out db/nomatching_[nom_souche].txt -outfmt 2

Ou, pour exécuter sur un cluster :

css

    qsub -cwd -V -N Blast.[nom_souche] -o qlogs -e qlogs -pe thread 20 -b y "conda activate blast-2.13.0 && blastn -db db/gbct_combined -query db/[nom_souche]_genomic.fna -out db/nomatching_[nom_souche].txt -outfmt 2 && conda deactivate"

Étape 6 : Extraction et Analyse des Séquences Uniques

    Lancer les scripts Python pour identifier les séquences uniques et celles ayant peu de correspondances

python unique_seq.py
python few_matches.py




## Étape 2: Création de la Base de Données

1. **Télécharger les Génomes RefSeq d'une Espèce Spécifique**
```
esearch -db assembly -query "[Nom_Espèce][Organism] AND latest_refseq[filter]" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F'/' '{print $0"/"$NF"_genomic.fna.gz"}' | sed 's|ftp://|https://|' | xargs -n 1 wget -P [Nom_Espèce]_genomes/
```

2. **Compter le Nombre de Génomes Téléchargés**
```
ls [Nom_Espèce]_genomes/*.fna.gz | wc -l
```

## Étape 3: Filtration et Préparation des Données

1. **Extraire et Analyser les Données du Génome**
```
gunzip -c gbct_[nom_souche]/[nom_fichier_génome].fna.gz | head -1
```

2. **Filtrer les Fichiers Génomiques**
```
{ find [Nom_Espèce]_genomes/ -name "*.fna.gz" -print0 | xargs -0 -I{} bash -c 'if gunzip -c "{}" | head -1 | grep -q "[motif_recherche]"; then echo "{}"; gunzip -c "{}" | head -1; fi'; } | tee >(grep -v '^>' | wc -l | xargs -I{} echo "Nombre de fichiers = {}")
```

## Étape 4: Construction de la Base de Données pour BLAST

1. **Créer le Répertoire de la Base de Données**
```
mkdir db
```

2. **Constituer la Base de Données**
```
for file in [Nom_Espèce]genomes/*.fna.gz; do gzip -dc "$file"; done > db/gbct[Nom_Espèce].fna
```

3. **Construire la Base de Données avec makeblastdb**
```
conda activate [env_blast] && makeblastdb -in db/gbct_[Nom_Espèce].fna -out db/gbct_[Nom_Espèce] -parse_seqids -dbtype nucl
```

## Étape 5: Réalisation du BLAST

1. **Exécuter le BLAST**
```
qsub -cwd -V -N Blast.[nom_souche] -o qlogs -e qlogs -pe thread 20 -b y "conda activate [env_blast] && blastn -db db/gbct_[Nom_Espèce] -query db/[nom_fichier_query].fna -out db/nomatching_[nom_souche].txt -outfmt 2 && conda deactivate"
```

2. **Analyse des Résultats**
- Lancer le script python <uniq_seq.py> pour l'analyse des séquences uniques.
- Conception des amorces/ primers de PCR (voir étape 6)


## Protocole 2 : Identification de primers spécifiques à l'échelle de l'espèce

Ce protocole détaille les étapes pour identifier des primers spécifiques à l'échelle de l'espèce bactérienne.

## Prérequis

- Système ou sous-système Linux (WSL) avec Python 3 installé.
- Environnements Conda avec BLAST+, Entrez-Direct, Prokka et Roary installés.

## Étape 1: Préparation de l'Environnement de Travail

1. **Créer un Répertoire de Travail et s'y déplacer**
```
mkdir -p [Nom_Espèce]primers/gbct[Nom_Espèce]
cd [Nom_Espèce]_primers
```

2. **Créer le Répertoire de la Base de Données**
```
mkdir db
```

## Étape 2: Téléchargement des Génomes de l'espèce d'intérêt

1. **Rechercher l'organisme "Muribaculum intestinale" dans la base de données RefSeq de NCBI et récupérer les liens FTP des génomes de référence**
esearch -db assembly -query "Muribaculum intestinale[Organism] AND latest_refseq[filter]" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F'/' '{print $0"/"$NF"_genomic.fna.gz"}' | sed 's|ftp://|https://|'


    Rechercher et Télécharger les Génomes RefSeq de l'Espèce Spécifique
    esearch -db assembly -query "[Nom_Espèce][Organism] AND latest_refseq[filter]" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F'/' '{print $0"/"$NF"_genomic.fna.gz"}' | sed 's|ftp://|https://|' | xargs -n 1 wget -P [Nom_Espèce]_genomes/

    Compter le Nombre de Génomes Téléchargés
    ls [Nom_Espèce]_genomes/*.fna.gz | wc -l

## Étape 3: Préparation des Données

    Afficher la Première Ligne de Chaque Fichier FASTA
    for file in [Nom_Espèce]_genomes/*.fna.gz; do echo -n "$file: "; zcat "$file" | head -1; done

    Identifier et Supprimer les Génomes Spécifiques (si nécessaire)
    for file in [Nom_Espèce]_genomes/*.fna.gz; do echo $file && gunzip -c $file | grep "[Motif_Recherche]"; done

## Étape 4: Annotation et Analyse des Génomes

    Annoter les Génomes avec Prokka
    for fna in [Nom_Espèce]genomes/*.fna; do qsub -cwd -V -N Prokka$(basename $fna .fna) -o qlogs -e qlogs -pe thread 20 -b y "conda activate prokka-x.x.x && prokka --outdir [Nom_Espèce]_annotations/$(basename $fna .fna) --prefix $(basename $fna .fna) --genus [Genus] --species [Species] --kingdom Bacteria --gcode 11 --cpus 20 $fna && conda deactivate"; done

    Utiliser Roary pour l'Analyse du Pangénome
    qsub -cwd -V -N Roary_[Nom_Espèce] -o qlogs.Roary -e qlogs.Roary -pe thread 40 -b y "conda activate roary-x.x.x && roary -p 40 -f Roary_[Nom_Espèce] -e --mafft -n -v -cd 100 [Nom_Espèce]_annotations//.gff && conda deactivate"

    Générer les Résultats et les Visualiser
    roary_plots.py --format png Roary_[Nom_Espèce]/core_gene_alignment.newick Roary_[Nom_Espèce]/gene_presence_absence.csv

    Extraire les Séquences du Core Genome
    cp Roary_[Nom_Espèce]/pan_genome_reference.fa db/[Nom_Espèce]_core_genome.fa

## Étape 5: Réalisation et Analyse du BLAST

    Réaliser le BLAST entre le Core Genome et la Base de Données
    blastn -db db/gbct_Bacteroidia -query db/[Nom_Espèce]core_genome.fa -out db/nomatching[Nom_Espèce]_core_genome.txt -outfmt 2

    Analyse des Résultats
        Lancer le script Python pour extraire les séquences uniques.
        Utiliser un outil de conception de primers PCR pour évaluer les primers potentiels.


## Étape 6 : Conception des amorces/ primers
- Examiner les fichiers de sortie pour identifier des primers potentiels à l'aide de l'outil "PCR Primer Design" (Eurofins) : https://eurofinsgenomics.eu/en/ecom/tools/pcr-primer-design/

Voici les paramètres et valeurs utilisés fréquemment (à adapter selon votre situation) : 

![image](https://github.com/MaximeNaour/GenoTools/blob/main/images/criteria.png?raw=true)
