# Protocole d'Identification de Primers Spécifiques

Ce protocole détaille les étapes pour identifier des primers spécifiques à l'échelle de la souche bactérienne (Protocole 1) et à l'échelle de l'espèce bactérienne (Protocole 2). 
## Prérequis

- Système ou sous-système Linux (WSL) avec Python 3 installé.
- Environnements Conda avec BLAST+ et Entrez-Direct installés.

## Protocole 1 : Identification de primers spécifiques d'une souche bactérienne

## Étape 1: Préparation de l'Environnement de Travail

1. **Créer un Répertoire de Travail**
```
mkdir -p [nom_souche]primers/gbct[nom_souche]
cd [nom_souche]_primers
```
2. **Télécharger le Génome de Référence**
```
wget [lien_ftp_génome_souche] -O gbct_[nom_souche]/[nom_fichier_génome].fna.gz
```

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
- Examiner les fichiers de sortie pour identifier des primers potentiels à l'aide de l'outil "PCR Primer Design" (Eurofins) : https://eurofinsgenomics.eu/en/ecom/tools/pcr-primer-design/

Voici les paramètres et valeurs utilisés fréquemment (à adapter selon votre situation) : 

<div align="center">
  ![unnamed](https://github.com/MaximeNaour/GenoTools/assets/55536880/b144fe81-5b5b-499a-b77a-add95e512af0)
</div>
