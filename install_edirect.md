# Installer Edirect (NCBI) sur Linux

#### Edirect est une suite d'outils permettant d'accéder aux bases de données du National Center for Biotechnology Information (NCBI) via la ligne de commande. Suivez les étapes ci-dessous pour installer Edirect sur votre système Linux.

### Étape 1 : Mise à jour du système

Avant de commencer, mettez à jour votre système Linux pour vous assurer que vous disposez des dernières versions de tous les paquets et dépendances nécessaires. Exécutez la commande suivante pour mettre à jour votre système :

```
sudo apt-get update && sudo apt-get upgrade
```

### Étape 2 : Installer les dépendances

Edirect nécessite quelques dépendances (notamment Perl) pour fonctionner correctement. Installez les dépendances suivantes en utilisant la commande ci-dessous :

```
sudo apt-get install -y perl cpanminus libxml-simple-perl libwww-perl libnet-perl libexpat1-dev libssl-dev libxml2-dev libjson-perl libdigest-md5-perl
```

### Etape 3 : Installer Anaconda3 dans votre environnement de travail pour utiliser conda

Se déplacer dans le répertoire HOME, créer un répertoire "pkg" pour contenir les différents outils et scripts nécessaire à votre travail et télécharger le script d'installation d'Anaconda3 (lien pour voir les différentes versions d'Anaconda : https://repo.anaconda.com/archive/)
```
cd ~ && mkdir pkg && cd pkg && wget https://repo.continuum.io/archive/Anaconda3-2023.07-1-Linux-x86_64.sh
```
Le tilde "~" permet d'accéder directement au répertoire HOME de l'utilisateur.  

Installation d'Anaconda3 : Se placer dans le répertoire "pkg" où se trouve l'exécutable .sh d'Anaconda3 téléchargé précédemment. 

```
cd ~ && cd pkg && bash Anaconda3-2023.07-1-Linux-x86_64.sh
```

Une fois l'installation lancée, différentes instructions vous seront demandées : 

  1. Appuyer sur la touche "ENTER" de votre clavier 
  2. Choisir le répertoire où sera installée Anaconda3. Pour ma part, j'ai changé la proposition faite par le script par ce chemin "/home/maxime-inrae/pkg/anaconda3" pour avoir mon environnement conda dans le répertoire "pkg"
  3. Initialisation de Anaconda3 en lançant la commande "conda init", veuillez écrire "yes" pour procéder à l'initialisation
  4. Réinitialisez l'invite de commande Linux

Pour vous assurez que votre environnement conda a bien été créé, lancez la commande suivante dans votre invite de commande. 
```
conda info --envs
```
Voux devriez obtenir ce type de résultat :  
```
# conda environments:
#
base                  *  /home/maxime-inrae/pkg/anaconda3
```

Pour activer votre environnement conda afin d'installer et utiliser des outils et scripts informatiques, lancez cette commande.
```
conda activate base
```
Pour information, la commande suivante permet de désactiver ou de sortir de l'environnement conda.
```
conda deactivate
```

Une fois votre environnement conda activé, vous pouvez lister tous les outils informatiques installés dans cet environnement conda. Pour ce faire, veuillez exécuter la commande suivante.
```
# packages in environment at /home/maxime-inrae/pkg/anaconda3:
#
# Name                    Version                   Build  Channel
_anaconda_depends         2023.07                 py311_0
...
```


### Étape 4 : Installer Edirect avec conda

Après avoir installé les dépendances nécessaires, vous pouvez maintenant procéder à l'installation d'Edirect. Exécutez la commande suivante pour télécharger et le script d'installation d'Edirect dans votre environnement conda :

```
conda activate base
conda create -y -n edirect -c conda-forge -c bioconda -c defaults -c defaults entrez-direct
```

Une fois l'installation faite, lancez la commande suivante pour s'assurer que l'outil EDIRECT a bien été installé dans votre environnement de travail conda.  
```
conda info --envs
```
Voux devriez obtenir ce type de résultat : 
```
# conda environments:
#
base                  *  /home/maxime-inrae/pkg/anaconda3
edirect                  /home/maxime-inrae/pkg/anaconda3/envs/edirect
```
Nous constatons bien que l'outil "EDIRECT" a bien été installé dans le dossier "envs/" du répertoire de travail Anaconda3.

### Étape 5 : Vérifier l'installation

Pour vérifier que l'outil "EDIRECT" a été correctement installé, exécutez la commande suivante :

```
esearch -version
```

Si l'installation a réussi, vous devriez voir le numéro de version d'Edirect s'afficher (ex: version 16.2 à la date de 31/07/2023)

Vous pouvez désormais commencer à utiliser Edirect pour interagir avec les bases de données du NCBI via la ligne de commande.
