///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

Marche à suivre pour Mepp sous Mac OS X ('10.5 Leopard'* ou '10.6 Snow Leopard') :
----------------------------------------------------------------------------------

1) installer Xcode:
http://developer.apple.com/technologies/xcode.html

2) installer MacPorts:
http://www.macports.org/install.php

3) mettre à jour MacPorts: 
sudo port -v selfupdate

4) installer les paquets suivants avec MacPorts: 
sudo port install cgal						(assez long, <= 1 heure)
sudo port install qt4-mac					(très long, environ 2-3 heures)
sudo port install libQGLViewer				(très rapide, quelques minutes)
sudo port install glew						(très très rapide, quelques secondes)

5) télécharger les sources de Mepp:
svn checkout https://nom-du-développeur@scm.gforge.liris.cnrs.fr/svnroot/mepp/trunk

6a) compiler Mepp avec CMake et Makefile:
se positionner dans trunk, puis,
 - mkdir build; cd build
 - pour une version Release: cmake .. -DCMAKE_C_FLAGS='-arch x86_64' puis make
 - pour une version Debug: cmake .. -DCMAKE_C_FLAGS='-arch x86_64' -DCMAKE_BUILD_TYPE=Debug (2 fois pour la prise en compte du mode Debug) puis make
note: vous pouvez également utiliser la version "graphique" de CMake: cmake-gui

ou

6b) compiler Mepp avec CMake et Qt Creator:
CMake (ou CMake-gui) est capable de générer (pour le moment) des Makefiles Unix ainsi que des projets Code::Blocks ou encore KDevelop,
mais malheureusement pas encore des projets Qt Creator (.pro).
Qt Creator sait par contre interprêter directement le CMakeLists.txt de cmake
pour générer les projets Qt Creator (.pro).

Donc, pour générer les projets Qt Creator (.pro), voici la marche à suivre au sein de Qt Creator:

Ouvrir un fichier et choisir en bas 'Fichier de projet CMake'
puis choisir le CMakeLists.txt de la racine du trunk puis cliquer sur 'suivant'
puis dans 'Arguments' mettre -DCMAKE_C_FLAGS="-arch x86_64" -DCMAKE_BUILD_TYPE=Debug pour générer les projets Qt Creator (.pro) en Debug
ou -DCMAKE_C_FLAGS="-arch x86_64" pour générer les projets Qt Creator (.pro) en Release.
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
"make install" permet de déployer un bundle portable sur un autre système Mac, le bundle généré se trouve dans le dossier mepp_deploy



* Sous Mac OS X 10.5 Leopard, gcc 4 contient un bug lors de la compilation optimisée d'une application avec la librairie CGAL,
il faut donc utiliser gcc 4.2 à la place:
cmake .. -DCMAKE_CXX_COMPILER=/usr/bin/g++-4.2 -DCMAKE_C_COMPILER=/usr/bin/gcc-4.2

De plus, si votre processeur ne supporte pas le 64 bits, il faut également passer -DCMAKE_C_FLAGS='-arch i386'
