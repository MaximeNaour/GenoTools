Voici un tutoriel étape par étape pour installer Ubuntu sur Windows à l'aide de Windows Subsystem for Linux (WSL) et Visual Studio Code.

Installation de Visual Studio Code et du Remote-WSL Extension

Étape 1: Installation de Visual Studio Code

Aller à https://code.visualstudio.com/
Téléchargez et installez la version appropriée pour votre système (Windows dans ce cas).
Une fois l'installation terminée, ouvrez Visual Studio Code.

Étape 2: Installation de Remote-WSL Extension

Dans Visual Studio Code, ouvrez le volet Extensions en cliquant sur le carré dans la barre d'activité sur le côté.
Recherchez "Remote - WSL".
Cliquez sur "Install" pour installer l'extension.

Associer une distribution Linux à Visual Studio Code

Etape 3: Installer une distribution Linux

Dans le volet gauche de VSC, cliquer sur le logo "Remote Explorer"
Sélectionner "WSL Targets" dans le menu déroulant puis cliquer sur "+" (Add a distro)
Choisir la distribution souhaitée parmis celles proposées dans le Microsoft Store Apps (exemple : Ubuntu)
Cliquer sur "Get in Store App" puis "Obtenir" et ensuite "Ouvrir"
Une invite de commandes devrait s'ouvrir vous demandant de renseigner un nouvel identifiant et un nouveau mot de passe. 

Utilisation de Visual Studio Code avec WSL

Étape 4: Ouvrir la distribution Linux dans Visual Studio Code

Dans le volet gauche de VSC, cliquer sur le logo "Remote Explorer"
Sélectionner "WSL Targets" dans le menu déroulant
Cliquer sur "Connect to WSL" (symbole fenêtre avec un "+")

Étape 5: Sélectionner son répertoire de travail
Dans le volet gauche de VSC, cliquer sur le logo "Explorer" (Ctrl+Shift+E) puis "Open Folder"
Pour travailler directement dans votre environnement Windows à l'aide du système Linux installé, cliquer sur "Show Local" et sélectionner le dossier/répertoire souhaité.

Félicitations, vous avez installé Ubuntu sur Windows à l'aide de WSL et vous pouvez maintenant y accéder via Visual Studio Code ! 
Vous pouvez maintenant exécuter des commandes Linux directement dans le terminal intégré de Visual Studio Code.

