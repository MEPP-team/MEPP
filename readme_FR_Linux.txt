///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010-2011
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

----------------------------------------------------------------
(A) Marche à suivre pour Mepp sous Ubuntu 9.10, 10.04 ou 10.10 :
----------------------------------------------------------------

1) installer Ubuntu 9.10, 10.04 ou 10.10
2) faire la mise à jour de la distribution et des paquets déjà installés:
sudo apt-get update
3) installer les paquets suivants:
sudo apt-get install subversion libcgal-dev qtcreator libqglviewer-qt4-dev g++ cmake-gui libglew-dev doxygen-gui graphviz libxerces-c-dev
sudo apt-get install mencoder ffmpeg unrar filezilla

-------------------------------
Suivant votre version d'Ubuntu:
-------------------------------
Sous Ubuntu 9.10 créer 2 liens:
cd /usr/lib
sudo ln -s libqglviewer-qt4.so libQGLViewer.so
cd /usr/include
sudo ln -s qglviewer-qt4 QGLViewer

ou

sous Ubuntu 10.04 créer 1 lien:
cd /usr/lib
sudo ln -s libqglviewer-qt4.so libQGLViewer.so

ou

sous Ubuntu 10.10, rien à faire
-------------------------------


----------------------------------------------
(B) Marche à suivre pour Mepp sous Fedora 14 : http://pkgs.org/
----------------------------------------------

1) installer Fedora 14
2) faire la mise à jour de la distribution et des paquets déjà installés
3) installer les paquets suivants:
sudo yum install subversion CGAL-devel qt-creator libQGLViewer-devel gcc-c++ make cmake-gui glew-devel doxygen graphviz xerces-c-dev


----------------------------------------------
----------------------------------------------


4) télécharger les sources de Mepp:
svn checkout https://nom-du-développeur@scm.gforge.liris.cnrs.fr/svnroot/mepp/trunk

5a) compiler Mepp avec CMake et Makefile:
se positionner dans trunk, puis,
 - mkdir build; cd build
 - pour une version Release: "cmake .." puis make
 - pour une version Debug: "cmake .. -DCMAKE_BUILD_TYPE=Debug" (2 fois pour la prise en compte du mode Debug) puis make
note: vous pouvez également utiliser la version "graphique" de CMake: cmake-gui

 - la documentation de Mepp (à venir...) se génère avec: make mepp_DOC
 - la documentation de votre composant au format Doxygen se génère avec: make component_nomducomposant_DOC (exemple: make component_CGAL_Example_DOC)

ou

5b) compiler Mepp avec CMake et Qt Creator:
CMake (ou CMake-gui) est capable de générer (pour le moment) des Makefiles Unix ainsi que des projets Code::Blocks ou encore KDevelop,
mais malheureusement pas encore des projets Qt Creator (.pro).
Qt Creator sait par contre interprêter directement le CMakeLists.txt de cmake
pour générer les projets Qt Creator (.pro).

Donc, pour générer les projets Qt Creator (.pro), voici la marche à suivre au sein de Qt Creator:

Ouvrir un fichier et choisir en bas 'Fichier de projet CMake'
puis choisir le CMakeLists.txt de la racine du trunk puis cliquer sur 'suivant'
puis dans 'Arguments' mettre -DCMAKE_BUILD_TYPE=Debug pour générer les projets Qt Creator (.pro) en Debug
ou rien pour générer les projets Qt Creator (.pro) en Release.
Cliquer ensuite sur 'Exécuter CMake', puis, après le déroulement de CMake, cliquer sur 'Terminer' et le projet (.pro) se charge dans l'IDE.

----------------------------------------------------------------------
Note au sujet de l'activation/désactivation des composants avec CMake:
----------------------------------------------------------------------
- activation des composants: passer -DBUILD_component_nomducomposant=ON à CMake (exemple: -DBUILD_component_Curvature=ON)
- désactivation des composants: passer -DBUILD_component_nomducomposant=OFF à CMake (exemple: -DBUILD_component_Curvature=OFF)
