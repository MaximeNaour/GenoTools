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

### Étape 3 : Installer Edirect

Après avoir installé les dépendances nécessaires, vous pouvez maintenant procéder à l'installation d'Edirect. Exécutez la commande suivante pour télécharger le script d'installation d'Edirect :

```
cd ~ && wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz
```

Maintenant, décompressez l'archive téléchargée :

```
tar -xzf edirect.tar.gz
```

Accédez au répertoire nouvellement créé et exécutez le script d'installation :

```
cd edirect && ./setup.sh
```

### Étape 4 : Ajouter Edirect à votre PATH

Le fichier .bashrc se trouve normalement dans le chemin suivant "/home/\<username>/". Voici un exemple d'une commande bash pour s'assurer que le fichier .bashrc se trouve bien à cet emplacement (Remplacez \<maxime-inrae> par votre identifiant indiqué dans le dossier /home) : 

```
ls -lah /home/maxime-inrae/
```

Pour pouvoir utiliser Edirect depuis n'importe quel répertoire, vous devez l'ajouter à votre variable d'environnement PATH. Pour ce faire, exécutez la commande suivante :

```
echo 'export PATH=$HOME/edirect:$PATH' >> ~/.bashrc
```

Rechargez le fichier .bashrc pour que les modifications prennent effet :

```
source ~/.bashrc
```

### Étape 5 : Vérifier l'installation

Pour vérifier si Edirect a été correctement installé, exécutez la commande suivante :

```
esearch -version
```

Si l'installation a réussi, vous devriez voir le numéro de version d'Edirect s'afficher.

Vous pouvez désormais commencer à utiliser Edirect pour interagir avec les bases de données du NCBI via la ligne de commande.
