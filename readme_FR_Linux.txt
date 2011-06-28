///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010-2011
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

---------------------------------------------------------------------------------
(A) Marche à suivre pour Mepp sous Ubuntu 9.10, 10.04, 10.10 ou 11.04 / Debian 6:
---------------------------------------------------------------------------------

1) installer Ubuntu 9.10, 10.04, 10.10 ou 11.04 / Debian 6
2) faire la mise à jour de la distribution et des paquets déjà installés:
sudo apt-get update
3) installer les paquets suivants:
sudo apt-get install subversion libcgal-dev qtcreator libqglviewer-qt4-dev g++ cmake-gui libglew-dev doxygen graphviz libxerces-c-dev
sudo apt-get install libavcodec-dev libavformat-dev libavdevice-dev libswscale-dev

----------------------------------------
Suivant votre version d'Ubuntu / Debian:
----------------------------------------
Sous Ubuntu 9.10 (The Karmic Koala - le koala karmique), créer 2 liens:
 cd /usr/lib
 sudo ln -s libqglviewer-qt4.so libQGLViewer.so
 cd /usr/include
 sudo ln -s qglviewer-qt4 QGLViewer

ou

sous Ubuntu 10.04 (The Lucid Lynx - le lynx lucide), créer 1 lien:
 cd /usr/lib
 sudo ln -s libqglviewer-qt4.so libQGLViewer.so

ou

sous Ubuntu 10.10 (The Maverick Meerkat - le suricate rebelle) / Debian 6, rien à faire

ou

sous Ubuntu 11.04 (The Natty Narwhal - le narval chic), créer 2 liens:
 cd /usr/lib
si Linux 32 bit
 sudo ln -s i386-linux-gnu/libX11.so libX11.so
 sudo ln -s i386-linux-gnu/libXext.so libXext.so
si Linux 64 bit
 sudo ln -s x86_64-linux-gnu/libX11.so libX11.so
 sudo ln -s x86_64-linux-gnu/libXext.so libXext.so
 
ou

sous Ubuntu 11.10 (The Oneiric Ocelot - l'ocelot onirique, testé avec 'pré-version' du 24 juin -> http://cdimage.ubuntu.com/daily-live/current/), créer 2 liens:
 cd /usr/lib
si Linux 32 bit
 sudo ln -s i386-linux-gnu/libGLU.so libGLU.so
 sudo ln -s i386-linux-gnu/libGL.so libGL.so
si Linux 64 bit
 sudo ln -s x86_64-linux-gnu/libGLU.so libGLU.so
 sudo ln -s x86_64-linux-gnu/libGL.so libGL.so
------------------------------------------


----------------------------------------------------
(B) Marche à suivre pour Mepp sous Fedora 14 ou 15 :
----------------------------------------------------

1) installer Fedora 14 ou 15
2) faire la mise à jour de la distribution et des paquets déjà installés (facultatif)
3) installer les paquets suivants:
sudo yum install subversion CGAL-devel qt-creator libQGLViewer-devel cmake-gui glew-devel doxygen graphviz xerces-c-devel


--------------------------------------------------
(C) Marche à suivre pour Mepp sous openSUSE 11.4 :
--------------------------------------------------

1) installer openSUSE 11.4
2) faire la mise à jour de la distribution et des paquets déjà installés (facultatif)
3a) ajouter le repository suivant à YaST2: http://download.opensuse.org/repositories/openSUSE:/11.4:/Contrib/standard/
3b) installer les paquets suivants:
sudo zypper install subversion CGAL-devel qt-creator cmake-gui glew-devel doxygen graphviz libxerces-c-devel

3c) télécharger libQGLViewer: http://www.libqglviewer.com/src/libQGLViewer-2.3.9.tar.gz
3d) compiler et installer libQGLViewer:
tar -xzf libQGLViewer-2.3.9.tar.gz
cd libQGLViewer-2.3.9/QGLViewer
qmake
make
sudo make install


--------------------------------------
(D) Pour toutes les versions de Linux:
--------------------------------------


4) télécharger les sources de Mepp:
svn checkout https://nom-du-développeur@scm.gforge.liris.cnrs.fr/svnroot/mepp (en prenant bien soin de renseigner "votre_username_gforge" pour les membres du LIRIS) ou
svn checkout http://scm.gforge.liris.cnrs.fr/public/mepp (pour les utilisateurs anonymes)

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
