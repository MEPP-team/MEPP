///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

Marche à suivre pour Mepp sous Mac OS X (10.5 ou 10.6) :
--------------------------------------------------------

1) installer Xcode:
http://developer.apple.com/technologies/xcode.html

2) installer MacPorts:
http://www.macports.org/install.php

3) mettre à jour MacPorts: 
sudo port -v selfupdate

4) installer les paquets suivants avec MacPorts: 
sudo port -v install cgal +universal
sudo port -v install qt4-mac +universal
sudo port -v install libQGLViewer +universal

5) télécharger les sources de Mepp:
svn checkout https://nom-du-développeur@scm.gforge.liris.cnrs.fr/svnroot/mepp/trunk

6a) compiler Mepp avec CMake et Makefile:
se positionner dans trunk, puis,
 - pour une version Release: cmake . -DCMAKE_C_FLAGS='-arch x86_64 -arch i386' puis make
 - pour une version Debug: cmake . -DCMAKE_C_FLAGS='-arch x86_64 -arch i386' -DCMAKE_BUILD_TYPE=Debug puis make
note: vous pouvez également utiliser la version "graphique" de CMake: cmake-gui

ou

6b) compiler Mepp avec CMake et Qt Creator:
CMake (ou CMake-gui) est capable de générer (pour le moment) sous Linux des Makefiles Unix ainsi que des projets Code::Blocks ou encore KDevelop,
mais malheureusement pas encore des projets Qt Creator (.pro).
Qt Creator sait par contre interprêter directement le CMakeLists.txt de cmake
pour générer les projets Qt Creator (.pro).

Donc, pour générer les projets Qt Creator (.pro), voici la marche à suivre au sein de Qt Creator:

Ouvrir un fichier et choisir en bas 'Fichier de projet CMake'
puis choisir le CMakeLists.txt de la racine du trunk puis cliquer sur 'suivant'
puis dans 'Arguments' mettre -DCMAKE_C_FLAGS='-arch x86_64 -arch i386' -DCMAKE_BUILD_TYPE=Debug pour générer les projets Qt Creator (.pro) en Debug
ou -DCMAKE_C_FLAGS='-arch x86_64 -arch i386' pour générer les projets Qt Creator (.pro) en Release.
Cliquer ensuite sur 'Exécuter CMake', puis, après le déroulement de CMake, cliquer sur 'Terminer' et le projet (.pro) se charge dans l'IDE.

----------------------------------------------------------------------
Note au sujet de l'activation/désactivation des composants avec CMake:
----------------------------------------------------------------------
- activation des composants: passer -DBUILD_component_nomducomposant=ON à CMake (exemple: -DBUILD_component_Curvature=ON)
- désactivation des composants: passer -DBUILD_component_nomducomposant=OFF à CMake (exemple: -DBUILD_component_Curvature=OFF)
