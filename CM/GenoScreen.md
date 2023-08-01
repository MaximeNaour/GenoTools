## *** Tutoriel spécial Catherine : A l'exploration des génomes ***

### Avant de partir à l'aventure, il est nécessaire de configurer notre espace de travail.  

#### Etape 1 : Configurer un environnement de travail avec une distribution Linux - Indispensable pour utiliser des commandes ou outils bioinformatiques.

Voir le fichier "wsl_VSC.md" pour travailler dans un environnement Linux/Ubuntu par le biai d'un Windows Subsystem Linux (WLS) et de Visual Studio Code (VSC).  

#### Etape 2 : Installer EDirect dans votre environnement de travail Linux  

Se référer au fichier "install_edirect.md" situé directement dans le répertoire "GenoTools".

La partie la plus fastidieuse est faite Catherine "*pfiou*"... Maintenant laisse-toi guider par ce tutoriel conçu rien que pour toi !

### Explorons les génomes !

#### Télécharger les génomes d'intérêt

esearch -db assembly -query "Adlercreutzia equolifaciens subsp. celatus[ORGN]" | efetch -format docsum > assemblies.xml  
grep -E '<SubmissionDate>|<FtpPath_GenBank>' assemblies.xml | awk 'NR%2{printf "%s ",$0;next;}1' | sed -E 's/\<SubmissionDate>([^<]+)<\/SubmissionDate> <FtpPath_GenBank>([^<]+)<\/FtpPath_GenBank>/\1,\2/' | sort -t',' -k1 -r | cut -d',' -f2  
esearch -db assembly -query "Adlercreutzia equolifaciens subsp. celatus[ORGN]" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F'/' '{print $0"/"$NF"_genomic.gbff.gz"}'  
esearch -db taxonomy -query "Adlercreutzia equolifaciens subsp. celatus" | efetch -format docsum | xtract -pattern DocumentSummary -element Id    
