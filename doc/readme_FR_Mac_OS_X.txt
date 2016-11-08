///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010-2011-2012-2013-2014-2015-2016
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

Marche à suivre pour Mepp sous Mac OS X ('10.5 Leopard'*, '10.6 Snow Leopard', '10.7 Lion', '10.8 Mountain Lion', '10.9 Mavericks', '10.10 Yosemite', '10.11 El Capitan' ou '10.12 Sierra') :
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

1a) A partir de Mac OS X '10.8 Mountain Lion', X11 ne fait plus partie intégrante du système, vous devez donc dans ce cas installer à la place XQuartz:
http://xquartz.macosforge.org/

---> avec une version récente de Homebrew (cf. point 3 ci-dessous), vous pouvez aussi installer XQuartz avec Homebrew Cask:
brew install Caskroom/cask/xquartz

1b) installer Xcode:
http://developer.apple.com/xcode/ ou via l'App Store

---> pour une version de Xcode >= 4.3 mais < 5.0, vous devez aussi installer les 'outils en ligne de commande' pour Xcode.
Vous pouvez les installer à l'intérieur de Xcode (Download preferences). Vous aurez aussi également besoin d'éxécuter cette commande:
sudo xcode-select -switch /Applications/Xcode.app/Contents/Developer

---> à partir de Mac OS X '10.9 Mavericks', l'installation des 'outils en ligne de commande' se fait désormais avec la commande suivante:
xcode-select --install

2) installer Java Developer Package for Mac OS X
https://connect.apple.com/
NOTE: ---> cela n'est pas utile à partir de Mac OS X '10.9 Mavericks'

///////////////////////////////////////////////////////////////////////////

3) utiliser Homebrew (3a) ou MacPorts (3b)
Je vous recommande fortement Homebrew, plus puissant et plus rapide lors de l'installation car disposant notamment parfois de "Bottle" (paquet binaire)

3a) installer Homebrew (http://brew.sh/):
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

---> mettre à jour Homebrew: 
brew update

---> facultatif:
brew tap homebrew/science

---> installer les paquets suivants avec Homebrew:
brew install cgal
brew install qt (sous '10.12 Sierra' : brew install cartr/qt4/qt)
brew install libqglviewer
brew install glew
brew install doxygen graphviz xerces-c
brew install ffmpeg
brew install assimp

brew install python

---> il vous sera peut-être utile de créer ce lien pour libqglviewer (attention ici, 2.6.3 est un exemple...)
sudo ln -s "/usr/local/Cellar/libqglviewer/2.6.3/lib/QGLViewer.framework" "/Library/Frameworks/QGLViewer.framework"

NOTE: ---> sous Mac OS X '10.9 Mavericks' uniquement, il vous faudra exécuter les 2 commandes suivantes:
sudo ln -s /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/System/Library/Frameworks/OpenGL.framework/Headers /System/Library/Frameworks/OpenGL.framework/Headers
sudo ln -s /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.9.sdk/System/Library/Frameworks/AGL.framework/Headers /System/Library/Frameworks/AGL.framework/Headers


ou


3b) installer MacPorts:
http://www.macports.org/install.php

---> mettre à jour MacPorts: 
sudo port -v selfupdate

---> installer les paquets suivants avec MacPorts:
sudo port install cgal												(assez long, <= 1 heure)
sudo port install qt4-mac											(très long, environ 2-3 heures)
sudo port install libQGLViewer										(très rapide, quelques minutes)
sudo port install glew												(très très rapide, quelques secondes)
sudo port install doxygen graphviz xercesc3
sudo port install ffmpeg											(très rapide, quelques minutes)

---> puis:
sudo port install python_select
sudo python_select pythonXX (ex.: python27)

///////////////////////////////////////////////////////////////////////////

4a) compiler Mepp avec CMake et Makefile:
se positionner dans le dossier MEPP.git (ou/or MEPP.svn)/MEPP, puis,
 - mkdir build; cd build
 - pour une version Release: "cmake .." puis make
 - pour une version Debug: "cmake .. -DCMAKE_BUILD_TYPE=Debug" (2 fois pour la prise en compte du mode Debug) puis make
note: vous pouvez également utiliser la version "graphique" de CMake: cmake-gui

 - la documentation de Mepp se génère avec: make mepp_DOC
 - la documentation de votre composant au format Doxygen se génère avec: make component_nomducomposant_DOC (exemple: make component_CGAL_Example_DOC)

ou

4b) compiler Mepp avec CMake et Qt Creator:
CMake (ou CMake-gui) est capable de générer (pour le moment) des Makefiles Unix ainsi que des projets Code::Blocks ou encore KDevelop,
mais malheureusement pas encore des projets Qt Creator (.pro).
Qt Creator sait par contre interprêter directement le CMakeLists.txt de cmake
pour générer les projets Qt Creator (.pro).

Donc, pour générer les projets Qt Creator (.pro), voici la marche à suivre au sein de Qt Creator:

Ouvrir un fichier et choisir en bas 'Fichier de projet CMake'
puis choisir le CMakeLists.txt du dossier MEPP.git (ou/or MEPP.svn)/MEPP (attention, pas le dossier MEPP.git (ou/or MEPP.svn)/MEPP/src !) puis cliquer sur 'suivant'
puis dans 'Arguments' mettre -DCMAKE_BUILD_TYPE=Debug pour générer les projets Qt Creator (.pro) en Debug
ou rien pour générer les projets Qt Creator (.pro) en Release.
Renseigner si besoin le chemin pour CMake: /opt/local/bin/cmake 
Cliquer ensuite sur 'Exécuter CMake', puis, après le déroulement de CMake, cliquer sur 'Terminer' et le projet (.pro) se charge dans l'IDE.

----------------------------------------------------------------------
Note au sujet de l'activation/désactivation des composants avec CMake:
----------------------------------------------------------------------
- activation des composants: passer -DBUILD_component_nomducomposant=ON à CMake (exemple: -DBUILD_component_Curvature=ON)
- désactivation des composants: passer -DBUILD_component_nomducomposant=OFF à CMake (exemple: -DBUILD_component_Curvature=OFF)

---------------------------------------------------------------
Note au sujet de la génération d'un bundle portable (mepp.app):
---------------------------------------------------------------
"sudo make install" permet de déployer un bundle portable sur un autre système Mac, le bundle généré se trouve dans le dossier mepp_deploy



* Sous 'Mac OS X 10.5 Leopard', gcc 4 contient un bug lors de la compilation optimisée d'une application avec la librairie CGAL,
il faut donc utiliser gcc 4.2 à la place:
cmake .. -DCMAKE_CXX_COMPILER=/usr/bin/g++-4.2 -DCMAKE_C_COMPILER=/usr/bin/gcc-4.2

De plus, si votre processeur ne supporte pas le 64 bits, il faut également passer -DCMAKE_C_FLAGS='-arch i386'
