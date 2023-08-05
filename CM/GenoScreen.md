<h2 align="center"><b>Tutoriel spécial Catherine : A l'exploration des génomes</b></h2>

<h3 align="center"><b>Etape 1 : Configuration de l'espace de travail</b></h3>

#### Avant de partir à l'aventure, il est nécessaire de configurer notre espace de travail.  

#### Partie 1 : Configurer un environnement de travail avec une distribution Linux - Indispensable pour utiliser des commandes ou outils bioinformatiques.

Voir le fichier "wsl_VSC.md" pour travailler dans un environnement Linux/Ubuntu par le biai d'un Windows Subsystem Linux (WLS) et de Visual Studio Code (VSC).  

#### Partie 2 : Installer EDirect dans votre environnement de travail Linux  

Se référer au fichier "install_edirect.md" situé directement dans le répertoire "GenoTools".

La partie la plus fastidieuse est faite Catherine "*pfiou*"... Maintenant laisse-toi guider par ce tutoriel conçu rien que pour toi !

<h3 align="center"><b>Etape 2 : Identification du génome d'intérêt avec Edirect : Précision à portée de main</b></h3>

La partie suivante est dédiée à l'utilisation de l'outil NCBI en ligne de commande "EDirect" pour collecter et télécharger les informations génomiques de la bactérie "Adlercreutzia equolifaciens subsp. celatus". Bien sûr, les lignes de commande suivantes peuvent être appliquées à tout autre organisme.

#### Avant de se lancer dans les commandes informatiques, il est nécessaire d'avoir installé sur son environnement de travail "conda" et "Edirect" (voir étape 1 - partie 2 si ces outils ne sont pas installés). 

##### Activer l'outil "Edirect" installé dans votre environnement conda. 

```
conda activate edirect
```

#### Partie 1 : Récupérer l'ID taxonomique pour "Adlercreutzia equolifaciens subsp. celatus" afin de télécharger son génome.

Le Taxonomy ID (taxid) est un identifiant unique du NCBI (National Center for Biotechnology Information) assigné à chaque entité taxonomique, garantissant l'unicité et évitant toute confusion entre les espèces à noms similaires. Il assure la stabilité en restant constant malgré les éventuels changements de nom d'un organisme au fil du temps. Le taxid permet également de déterminer les relations taxonomiques entre différents organismes et facilite la recherche de données dans les bases de données biologiques, ce qui s'avère plus précis et efficace que l'utilisation du nom de l'organisme.  

```
esearch -db taxonomy -query "Adlercreutzia equolifaciens subsp. celatus[ORGN]" | efetch -format docsum | xtract -pattern DocumentSummary -element Id    
```
- 'esearch' cherche dans la base de données de taxonomie pour "Adlercreutzia equolifaciens subsp. celatus".
- 'efetch' récupère le résumé du document
- 'xtract' est utilisé pour obtenir l'ID de taxonomie pour cet organisme

Résultat 2 taxids : 
1. 1121021 --> Adlercreutzia equolifaciens subsp. celatus DSM 18785 (souche)
2. 394340 --> Adlercreutzia equolifaciens subsp. celatus (sous-espèce)

#### Partie 2 : Télécharger les informations sur les assemblages de génomes pour "Adlercreutzia equolifaciens subsp. celatus"

```
esearch -db assembly -query "Adlercreutzia equolifaciens subsp. celatus[ORGN]" | efetch -format docsum > assemblies.xml  
```
ou à l'aide du taxid précédemment identifié :
```
esearch -db assembly -query "txid394340[Organism:exp]" | efetch -format docsum > assemblies.xml
```

- 'esearch' cherche dans la base de données "assembly" pour "Adlercreutzia equolifaciens subsp. celatus[ORGN]".
- 'efetch' télécharge les informations sous forme de résumé "docsum" et redirige le résultat dans le fichier "assemblies.xml".

#### Partie 3 : Extraire et organiser les dates de soumission et les chemins FTP à partir des informations du génome

```
grep -E '<Organism>|<SubmissionDate>|<FtpPath_GenBank>' assemblies.xml | awk 'NR%3{printf "%s ",$0;next;}1' | sed -E 's/<SpeciesName>([^<]+)<\/SpeciesName> <SubmissionDate>([^<]+)<\/SubmissionDate> <FtpPath_GenBank>([^<]+)<\/FtpPath_GenBank>/\1,\2,\3/' | sort -t',' -k2 -r | cut -d',' -f1,3   
```
- 'grep' est utilisé pour filtrer les noms des organismes, les dates de soumission et les chemins FTP à partir du fichier "assemblies.xml".
- 'awk' combine ces deux lignes en une seule.
- 'sed' est utilisé pour supprimer les balises XML et séparer les champs par une virgule
- 'sort' trie le résultat par date de soumission en ordre descendant
- 'cut' sélectionne la deuxième colonne, qui est le chemin FTP

Résultat : 

    <Organism>Adlercreutzia equolifaciens subsp. celatus DSM 18785 (actinobacteria)</Organism>     <SubmissionDate>2022/07/06 00:00</SubmissionDate>     <FtpPath_GenBank>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/024/171/685/GCA_024171685.1_ASM2417168v1</FtpPath_GenBank>
    <Organism>Adlercreutzia equolifaciens subsp. celatus DSM 18785 (actinobacteria)</Organism>     <SubmissionDate>2018/11/13 00:00</SubmissionDate>     <FtpPath_GenBank>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/726/015/GCA_003726015.1_ASM372601v1</FtpPath_GenBank>
    <Organism>Adlercreutzia equolifaciens subsp. celatus (actinobacteria)</Organism>     <SubmissionDate>2022/07/30 00:00</SubmissionDate>     <FtpPath_GenBank>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/934/882/975/GCA_934882975.1_MTG246_bin.28.fa</FtpPath_GenBank>
    <Organism>Adlercreutzia equolifaciens subsp. celatus (actinobacteria)</Organism>     <SubmissionDate>2022/07/30 00:00</SubmissionDate>     <FtpPath_GenBank>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/934/879/575/GCA_934879575.1_MTG247_bin.17.fa</FtpPath_GenBank>
    <Organism>Adlercreutzia equolifaciens subsp. celatus (actinobacteria)</Organism>     <SubmissionDate>2022/07/30 00:00</SubmissionDate>     <FtpPath_GenBank>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/934/878/345/GCA_934878345.1_MTG248_bin.40.fa</FtpPath_GenBank>
    <Organism>Adlercreutzia equolifaciens subsp. celatus (actinobacteria)</Organism>     <SubmissionDate>2022/07/30 00:00</SubmissionDate>     <FtpPath_GenBank>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/934/854/485/GCA_934854485.1_MTG235_bin.30.fa</FtpPath_GenBank>
    <Organism>Adlercreutzia equolifaciens subsp. celatus (actinobacteria)</Organism>     <SubmissionDate>2022/07/30 00:00</SubmissionDate>     <FtpPath_GenBank>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/934/851/295/GCA_934851295.1_MTG233_bin.41.fa</FtpPath_GenBank>
    <Organism>Adlercreutzia equolifaciens subsp. celatus (actinobacteria)</Organism>     <SubmissionDate>2022/07/30 00:00</SubmissionDate>     <FtpPath_GenBank>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/934/847/765/GCA_934847765.1_MTG234_bin.67.fa</FtpPath_GenBank>
    <Organism>Adlercreutzia equolifaciens subsp. celatus (actinobacteria)</Organism>     <SubmissionDate>2021/02/05 00:00</SubmissionDate>     <FtpPath_GenBank>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/016/865/385/GCA_016865385.1_ASM1686538v1</FtpPath_GenBank>
    <Organism>Adlercreutzia equolifaciens subsp. celatus (actinobacteria)</Organism>     <SubmissionDate>2018/08/26 00:00</SubmissionDate>     <FtpPath_GenBank>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/428/485/GCA_003428485.1_ASM342848v1</FtpPath_GenBank>
    <Organism>Adlercreutzia equolifaciens subsp. celatus (actinobacteria)</Organism>     <SubmissionDate>2018/07/25 00:00</SubmissionDate>     <FtpPath_GenBank>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/340/325/GCA_003340325.1_ASM334032v1</FtpPath_GenBank>
    <Organism>Adlercreutzia equolifaciens subsp. celatus (actinobacteria)</Organism>     <SubmissionDate>2018/07/25 00:00</SubmissionDate>     <FtpPath_GenBank>ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/340/305/GCA_003340305.1_ASM334030v1</FtpPath_GenBank>

Les chemins FTP pour accéder aux différents génomes assemblés pour l'espèce d'intérêt sont triés par ordre chronologique du plus récent au plus ancien. Cette étape nous permet d'identifier le génome d'intérêt, par exemple, celui-ci "GCA_024171685.1_ASM2417168v1" correspondant au plus récent génome assemblé (premier lien). Une fois le nom du génome d'intérêt trouvé, passons à l'étape suivant pour identifier le chemin 'ftp' puis 'https' du fichier '.gbff' d'intérêt.

Les fichiers .gbff sont des fichiers, de format GenBank Flat File, utilisés pour stocker des données de séquence génomique annotée. 

#### Partie 4 : Récupérer le chemin FTP et HTTPS pour le génome de référence

```
esearch -db assembly -query "Adlercreutzia equolifaciens subsp. celatus[ORGN]" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F'/' '{print $0"/"$NF"_genomic.gbff.gz"}' | sed 's|ftp://|https://|'    
```
ou à l'aide du taxid d'intérêt :
```
esearch -db assembly -query "txid394340[Organism:exp]" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F'/' '{print $0"/"$NF"_genomic.gbff.gz"}' | sed 's|ftp://|https://|'     
```

- 'xtract' est utilisé pour extraire le chemin FTP pour le génome de référence
- 'awk' ajoute le nom du fichier génomique à la fin du chemin FTP
- sed 's|ftp://|https://|' remplace 'ftp://' par 'https://'

Résultat : 

    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/934/879/575/GCF_934879575.1_MTG247_bin.17.fa/GCF_934879575.1_MTG247_bin.17.fa_genomic.gbff.gz
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/024/171/685/GCF_024171685.1_ASM2417168v1/GCF_024171685.1_ASM2417168v1_genomic.gbff.gz
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/865/385/GCF_016865385.1_ASM1686538v1/GCF_016865385.1_ASM1686538v1_genomic.gbff.gz
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/726/015/GCF_003726015.1_ASM372601v1/GCF_003726015.1_ASM372601v1_genomic.gbff.gz
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/428/485/GCF_003428485.1_ASM342848v1/GCF_003428485.1_ASM342848v1_genomic.gbff.gz
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/340/325/GCF_003340325.1_ASM334032v1/GCF_003340325.1_ASM334032v1_genomic.gbff.gz
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/340/305/GCF_003340305.1_ASM334030v1/GCF_003340305.1_ASM334030v1_genomic.gbff.gz

D'après le nom du génome d'intérêt, identifié précédemment, nous allons télécharger le fichier suivant : 

    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/024/171/685/GCF_024171685.1_ASM2417168v1/GCF_024171685.1_ASM2417168v1_genomic.gbff.gz

#### Partie 5 : Télécharger et décompresser le fichier d'annotation génomique GBFF

##### Télécharger le fichier dans le dossier spécifique à l'espèce d'intérêt
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/024/171/685/GCF_024171685.1_ASM2417168v1/GCF_024171685.1_ASM2417168v1_genomic.gbff.gz -P Ae_celatus
```

##### Se déplacer dans le dossier d'intérêt
```
cd Ae_celatus
```

##### Décompresser le fichier
```
gunzip GCF_024171685.1_ASM2417168v1_genomic.gbff.gz
```

##### Voir l'en-tête du fichier décompressé pour savoir si le fichier téléchargé est correct
```
head Ae_celatus/GCF_024171685.1_ASM2417168v1_genomic.gbff
```
Résultat : 
```
LOCUS       NZ_JAMTCE010000001    261222 bp    DNA     linear   CON 16-JUL-2023
DEFINITION  Adlercreutzia equolifaciens subsp. celatus DSM 18785
            BR46DRAFT_scaffold00001.1, whole genome shotgun sequence.
ACCESSION   NZ_JAMTCE010000001 NZ_JAMTCE010000000
VERSION     NZ_JAMTCE010000001.1
DBLINK      BioProject: PRJNA224116
            BioSample: SAMN02745888
            Assembly: GCF_024171685.1
KEYWORDS    WGS; RefSeq.
SOURCE      Adlercreutzia equolifaciens subsp. celatus DSM 18785
```

Si les informations renseignées te semblent correctes, nous pouvons passer à la prochaine étape du tutoriel qui vise à décortiquer les informations génomiques présentes dans le fichier 'GBFF'. Si les informations présentent dans ce fichier ne te convient pas, télécharge d'autres liens 'https' (voir 'partie 4' pour choisir un nouveau lien et 'partie 5' pour le télécharger et décompresser).

<h3 align="center"><b>### Etape 3 : Explorons le génome !</b></h3>

