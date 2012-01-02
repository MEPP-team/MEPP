///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010-2011-2012
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

Marche à suivre pour Mepp sous Mac OS X ('10.5 Leopard'*, '10.6 Snow Leopard' ou '10.7 Lion') :
-----------------------------------------------------------------------------------------------

1) installer Xcode:
http://developer.apple.com/xcode/

2) installer Java Developer Package for Mac OS X
https://connect.apple.com/

///////////////////////////////////////////////////////////////////////////

3) utiliser Homebrew (3a) ou MacPorts (3b)
Je vous recommande fortement Homebrew, plus puissant et plus rapide lors de l'installation car disposant notamment parfois de "Bottle" (paquet binaire)

3a) installer Homebrew (http://mxcl.github.com/homebrew/):
/usr/bin/ruby -e "$(curl -fsSL https://raw.github.com/gist/323731)"

---> !!! sous Mac OS X 10.7.x Lion avec Xcode 4.2, mettre à jour la formule Ruby de gmp (dépendance de cgal) !!! :
remplacer donc /usr/local/Library/Formula/gmp.rb par ce fichier: http://liris.cnrs.fr/mepp/download/gmp.rb
(ce point n'est pas nécessaire sous Mac OS X 10.7.x Lion avec Xcode 4.1)

---> installer les paquets suivants avec Homebrew:
brew install cgal
brew install qt
brew install glew
brew install doxygen graphviz xerces-c
brew install ffmpeg

---> télécharger libQGLViewer: http://www.libqglviewer.com/src/libQGLViewer-2.3.10.tar.gz
---> puis compiler et installer libQGLViewer:
tar -xzf libQGLViewer-2.3.10.tar.gz
cd libQGLViewer-2.3.10/QGLViewer
---> !!! supprimer la ligne 166 du fichier VRenderInterface.Qt4.ui !!!
qmake -spec macx-g++
sudo make install

ou

3b) installer MacPorts:
http://www.macports.org/install.php

---> mettre à jour MacPorts: 
sudo port -v selfupdate

---> installer les paquets suivants avec MacPorts:
sudo port install cgal				(assez long, <= 1 heure)
sudo port install qt4-mac			(très long, environ 2-3 heures)
sudo port install libQGLViewer		(très rapide, quelques minutes)
sudo port install glew				(très très rapide, quelques secondes)
sudo port install doxygen graphviz xercesc3
sudo port install ffmpeg			(très rapide, quelques minutes)

---> puis:
sudo port install python_select
sudo python_select pythonXX (ex.: python27)

///////////////////////////////////////////////////////////////////////////

4a) compiler Mepp avec CMake et Makefile:
se positionner dans le dossier MEPP.git/MEPP, puis,
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
puis choisir le CMakeLists.txt du dossier MEPP.git/MEPP (attention, pas le dossier MEPP.git/MEPP/src !) puis cliquer sur 'suivant'
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
