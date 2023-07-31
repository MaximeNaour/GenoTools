## Voici un tutoriel étape par étape pour installer une distribution Linux sur Windows à l'aide de Windows Subsystem for Linux (WSL) et Visual Studio Code.

### Installation de Visual Studio Code et du Remote-WSL Extension

#### Étape 1: Installation de Visual Studio Code
Aller à https://code.visualstudio.com/  
Téléchargez et installez la version appropriée pour votre système (Windows dans ce cas).  
Une fois l'installation terminée, ouvrez Visual Studio Code.  

#### Étape 2: Installation de WSL Extension
Dans Visual Studio Code, ouvrez le volet Extensions en cliquant sur le carré dans la barre d'activité sur le côté gauche. ![image](https://github.com/MaximeNaour/GenoTools/assets/55536880/21008143-d3d8-422c-9502-c9b32b70aa81)  
Recherchez "WSL".  
Cliquez sur "Install" pour installer l'extension.  

#### Étape 3: Activation de l'extension WSL
Pour activer "WSL" sur vore ordinateur Windows, veuillez suivre le lien suivant : https://pixiscreen.fr/ressources-dev/activer-le-sous-systeme-windows-pour-linux/  
##### Ne pas oublier de redémarrer l'ordinateur suite à l'activation. 

### Associer une distribution Linux à Visual Studio Code

#### Etape 4: Installer une distribution Linux
Dans le volet gauche de VSC, cliquer sur le logo "Remote Explorer" ![image](https://github.com/MaximeNaour/GenoTools/assets/55536880/43e2d994-4c58-492a-a068-2f76d72dc9e4)  
Sélectionner "WSL Targets" dans le menu déroulant puis cliquer sur "+" (Add a distro)  
Choisir la distribution souhaitée parmis celles proposées dans le Microsoft Store Apps (exemple : Ubuntu)  
Cliquer sur "Get in Store App" puis "Obtenir" et ensuite "Ouvrir"  
Une invite de commandes devrait s'ouvrir vous demandant de renseigner un nouvel identifiant et un nouveau mot de passe.   

### Utilisation de Visual Studio Code avec WSL

#### Étape 5: Ouvrir la distribution Linux dans Visual Studio Code
Dans le volet gauche de VSC, cliquer sur le logo "Remote Explorer" ![image](https://github.com/MaximeNaour/GenoTools/assets/55536880/43e2d994-4c58-492a-a068-2f76d72dc9e4)  
Sélectionner "WSL Targets" dans le menu déroulant  
Cliquer sur "Connect to WSL" (symbole fenêtre avec un "+")  

#### Étape 6: Sélectionner son répertoire de travail
Dans le volet gauche de VSC, cliquer sur le logo "Explorer" ![image](https://github.com/MaximeNaour/GenoTools/assets/55536880/1ee1c637-29b3-4f5f-8fba-9a5a39265bfd) (Ctrl+Shift+E) puis "Open Folder"  
Pour travailler directement dans votre environnement Windows à l'aide du système Linux installé, cliquer sur "Show Local" et sélectionner le dossier/répertoire souhaité.  

Félicitations, vous avez installé une distribution Linux sur Windows à l'aide de WSL et vous pouvez y accéder via Visual Studio Code!  
Vous pouvez maintenant exécuter des commandes Linux directement dans le terminal intégré de Visual Studio Code (Ctrl + ù).  
Pour changer de distribution système, cliquer sur la flèche à droite du "+" et sélectionner la distribution souhaitée.  

