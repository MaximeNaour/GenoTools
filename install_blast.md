# Tutoriel d'installation de BLAST+ sur WSL (Windows Subsystem for Linux)
### BLAST+ est un ensemble d'outils de ligne de commande pour effectuer des recherches BLAST (Basic Local Alignment Search Tool) sur des bases de données biologiques. 
## Ce tutoriel vous guidera pour installer BLAST+ sur le Windows Subsystem for Linux (WSL).

### Prérequis
Avoir installé le Windows Subsystem for Linux (WSL) avec une distribution Linux (par exemple, Ubuntu).  

### Étapes d'installation
Avant d'installer BLAST+, il est recommandé de mettre à jour les packages de votre système. 
```
sudo apt update && sudo apt upgrade
```

#### Téléchargement de BLAST+
Se déplacer dans un répertoire où vous stocker vos outils bioinformatiques  
```
mkdir -p $HOME/tools
cd $HOME/tools
```
Télécharger la dernière version de BLAST+ pour Linux  
```
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz 
```
Remplacez 2.15.0+ par la version actuelle que vous trouvez sur le site FTP (https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/).  

Extraction de l'archive  
```
tar -zxvpf ncbi-blast-2.15.0+-x64-linux.tar.gz
```
Cela créera un répertoire ncbi-blast-2.15.0+ dans votre répertoire courant.  

Ajout de BLAST+ à votre PATH  
Pour utiliser les outils BLAST+ dans n'importe quel répertoire, vous devez ajouter le répertoire bin à votre variable d'environnement PATH en modifiant le .bashrc  
```
echo 'export PATH=$PATH:$HOME/tools/ncbi-blast-2.15.0+/bin' >> ~/.bashrc
```
Appliquez les modifications  
```
source ~/.bashrc
```
Vérification de l'installation  
```
blastn -version
```

Vous devriez voir la version de BLAST+ affichée, ce qui confirme que l'installation a réussi.

### Configuration supplémentaire  
#### Il est possible de configurer un dossier pour stocker les bases de données BLAST

Créez un répertoire pour stocker vos bases de données BLAST  
```
mkdir $HOME/blastdb
```
Ajoutez la variable d'environnement BLASTDB à votre .bashrc  
```
echo 'export BLASTDB=$HOME/blastdb' >> ~/.bashrc
```
Appliquez les modifications  
```
source ~/.bashrc
```
Téléchargement des bases de données  
Utilisez le script update_blastdb.pl ($HOME/tools/ncbi-blast-2.15.0+/bin) pour télécharger les bases de données préformatées depuis NCBI  
Lister les bases de données BLAST disponible (veuillez attendre plusieurs secondes...)  
```
update_blastdb.pl --showall [*]
```
Exemple de commande : Télécharger et décompresser la base de données 16S_ribosomal_RNA dans le répertoire blastdb  
```
perl $HOME/tools/ncbi-blast-2.15.0+/bin/update_blastdb.pl --passive --decompress 16S_ribosomal_RNA
```
