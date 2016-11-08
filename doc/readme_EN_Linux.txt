///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010-2011-2012-2013-2014-2015-2016
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

------------------------------------------------------------------------------------------------------------------------------------------------------
(A) How-to for Mepp under Ubuntu 9.10, 10.04, 10.10, 11.04, 11.10, 12.04, 12.10, 13.04, 13.10, 14.04, 14.10, 15.04, 15.10, 16.04 or 16.10 / Debian 6 :
------------------------------------------------------------------------------------------------------------------------------------------------------

1) install Ubuntu 9.10, 10.04, 10.10, 11.04, 11.10, 12.04, 12.10, 13.04, 13.10, 14.04, 14.10, 15.04, 15.10, 16.04 or 16.10 / Debian 6
2) update distribution and packages:
sudo apt-get update
3) install new packages:
sudo apt-get install subversion libcgal-dev qtcreator libqglviewer-qt4-dev g++ cmake-gui libglew-dev doxygen graphviz libxerces-c-dev libassimp-dev
sudo apt-get install libavcodec-dev libavformat-dev libavdevice-dev libswscale-dev

NOTE1: ---> from Ubuntu 13.10 libqglviewer-qt4-dev package is renamed to libqglviewer-dev
NOTE2: ---> from Ubuntu 15.10 libqglviewer-qt4-dev package is renamed to libqglviewer-dev-qt4

NOTE3: ---> under Ubuntu 16.04 (no more problem under Ubuntu 16.10...) the libcgal-qt5-dev package is necessary (not normal, problem with the CGAL package !)

-------------------------------------------
Depending on your Ubuntu / Debian version :
-------------------------------------------
under Ubuntu 9.10 (The Karmic Koala), create 2 links:
 cd /usr/lib
 sudo ln -s libqglviewer-qt4.so libQGLViewer.so
 cd /usr/include
 sudo ln -s qglviewer-qt4 QGLViewer

or

under Ubuntu 10.04 (The Lucid Lynx), create 1 link:
 cd /usr/lib
 sudo ln -s libqglviewer-qt4.so libQGLViewer.so

or

under Ubuntu 10.10 (The Maverick Meerkat) / Debian 6, nothing to do

or

under Ubuntu 11.04 (The Natty Narwhal), create 2 links:
 cd /usr/lib
if Linux 32 bit
 sudo ln -s i386-linux-gnu/libX11.so libX11.so
 sudo ln -s i386-linux-gnu/libXext.so libXext.so
if Linux 64 bit
 sudo ln -s x86_64-linux-gnu/libX11.so libX11.so
 sudo ln -s x86_64-linux-gnu/libXext.so libXext.so
 
or

under Ubuntu 11.10* (The Oneiric Ocelot), create 5 links:
 cd /usr/lib
if Linux 32 bit
 sudo ln -s i386-linux-gnu/libQtOpenGL.so libQtOpenGL.so
 sudo ln -s i386-linux-gnu/libQtGui.so libQtGui.so
 sudo ln -s i386-linux-gnu/libQtCore.so libQtCore.so
 sudo ln -s i386-linux-gnu/libGLU.so libGLU.so
 sudo ln -s i386-linux-gnu/libGL.so libGL.so
if Linux 64 bit
 sudo ln -s x86_64-linux-gnu/libQtOpenGL.so libQtOpenGL.so
 sudo ln -s x86_64-linux-gnu/libQtGui.so libQtGui.so
 sudo ln -s x86_64-linux-gnu/libQtCore.so libQtCore.so
 sudo ln -s x86_64-linux-gnu/libGLU.so libGLU.so
 sudo ln -s x86_64-linux-gnu/libGL.so libGL.so

or

under Ubuntu 12.04* (The Precise Pangolin) and Ubuntu 12.10* (The Quantal Quetzal), nothing to do

or

under Ubuntu 13.04* (The Raring Ringtail), activate qt4 because by default qt5 is now active:
sudo apt-get install qt4-default

or

under Ubuntu 13.10* (The Saucy Salamander), activate qt4 because by default qt5 is now active:
sudo apt-get install qt4-default
and create 2 links:
 cd /usr/lib
if Linux 32 bit
 sudo ln -s i386-linux-gnu/libboost_thread.so libboost_thread.so
 sudo ln -s i386-linux-gnu/libboost_system.so libboost_system.so
if Linux 64 bit
 sudo ln -s x86_64-linux-gnu/libboost_thread.so libboost_thread.so
 sudo ln -s x86_64-linux-gnu/libboost_system.so libboost_system.so

or

under Ubuntu 14.04* (The Trusty Tahr), activate qt4 because by default qt5 is now active:
sudo apt-get install qt4-default

or

under Ubuntu 14.10* (The Utopic Unicorn), activate qt4 because by default qt5 is now active:
sudo apt-get install qt4-default

or

under Ubuntu 15.04* (The Vivid Vervet), activate qt4 because by default qt5 is now active:
sudo apt-get install qt4-default

or

under Ubuntu 15.10* (The Wily Werewolf), activate qt4 because by default qt5 is now active:
sudo apt-get install qt4-default

or

under Ubuntu 16.04* (The Xenial Xerus), activate qt4 because by default qt5 is now active:
sudo apt-get install qt4-default

or

under Ubuntu 16.04* (The Yakkety Yak), activate qt4 because by default qt5 is now active:
sudo apt-get install qt4-default


* to reinstall the standard GNOME desktop : sudo apt-get install gnome-panel
  to reinstall the Synaptic application :   sudo apt-get install synaptic


-----------------------------------------------
(B) How-to for Mepp under Fedora 14, 15 or 16 :
-----------------------------------------------

1) install Fedora 14, 15 or 16
2) update distribution and packages already installed (optional)
3) install new packages:
su (enter the password of the current user)
yum install subversion CGAL-devel qt-creator libQGLViewer-devel cmake-gui glew-devel doxygen graphviz xerces-c-devel


-----------------------------------------
(C) How-to for Mepp under openSUSE 11.4 :
-----------------------------------------

1) install openSUSE 11.4
2) update distribution and packages already installed (optional)
3a) add this repository to YaST2: http://download.opensuse.org/repositories/openSUSE:/11.4:/Contrib/standard/
3b) install new packages:
sudo zypper install subversion CGAL-devel qt-creator cmake-gui glew-devel doxygen graphviz libxerces-c-devel

3c) download libQGLViewer: http://www.libqglviewer.com/src/libQGLViewer-2.3.9.tar.gz
3d) compile and install libQGLViewer:
tar -xzf libQGLViewer-2.3.9.tar.gz
cd libQGLViewer-2.3.9/QGLViewer
qmake
make
sudo make install


----------------------------
(D) For all Linux versions :
----------------------------

4a) compile Mepp with CMake and Makefile:
go into MEPP.git (ou/or MEPP.svn)/MEPP folder, then,
 - mkdir build; cd build
 - for a Release version: "cmake .." then make
 - for a Debug version: "cmake .. -DCMAKE_BUILD_TYPE=Debug" (2 times for taking into account of debug mode) then make
note: you can also use the "graphics" version of CMake: cmake-gui

 - Mepp documentation can be generated with: make mepp_DOC
 - your component documentation in Doxygen format can be generated with: make component_componentname_DOC (example: make component_CGAL_Example_DOC)

or

4b) compile Mepp with CMake and Qt Creator:
CMake (or CMake-gui) is able to generate (right now) Unix Makefiles or Code::Blocks and KDevelop projects, but unfortunately not Qt Creator projects (.pro).
But fortunately Qt Creator can import directly CMakeLists.txt to generate Qt Creator projects (.pro).

So, in order to generate Qt Creator projects (.pro), here's the how-to under Qt Creator ("Applications > Programmation > Qt Creator" menu):

Open a file and choose at the bottom 'CMake project file'
then choose CMakeLists.txt in MEPP.git (ou/or MEPP.svn)/MEPP folder (warning, not MEPP.git (ou/or MEPP.svn)/MEPP/src folder !) then click on 'next'
then in 'Arguments' field put -DCMAKE_BUILD_TYPE=Debug for generate Qt Creator projects (.pro) in Debug mode
or nothing for generate Qt Creator projects (.pro) in Release mode.
Click then on 'Run CMake', then, after the execution of CMake, click on 'Terminate' and the project (.pro) loads into the IDE.

-------------------------------------------------------
Note on components activation/deactivation with CMake :
-------------------------------------------------------
- components activation: add -DBUILD_component_componentname=ON to CMake (example: -DBUILD_component_Curvature=ON)
- components deactivation: add -DBUILD_component_componentname=OFF to CMake (example: -DBUILD_component_Curvature=OFF)
