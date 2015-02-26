///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010-2011-2012-2013-2014-2015
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

How-to for Mepp (64 bits) under Windows with Visual Studio Express 2010 :
-------------------------------------------------------------------------


Note: you must have a 64-bit machine (Vista 64-bit or Seven 64-bit)

1a) download and install DAEMON Tools Lite in order to install the following 'iso' images:
http://download.gforge.liris.cnrs.fr/meppbin/windows/utils/ (DTLitexxxx-yyyy.exe)

1b) download and install Visual Studio Express 2010:
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/VS2010ExpressFRA.iso (728 Mo)

1c) download and install Windows SDK 7.1 (ie MSDK for Windows 7 and .NET Framework 4):
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/GRMSDKX_EN_DVD_7_and_4.0_(7.1).iso (571 Mo)
Note: if problem installing "MSDK for Windows 7 and .NET Framework 4" you can use web installer: http://www.microsoft.com/en-us/download/details.aspx?id=8279

1d) download and install Visual Studio Express 2010 SP1 (IMPORTANT):
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/sp1/VS2010SP1dvd1.iso (1.5 Go)

1e) download and install Visual C++ 2010 SP1 Compiler Update for the Windows SDK 7.1 (IMPORTANT):
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/sp1/VC-Compiler-KB2519277.exe (121 Mo)

2a) download MEPP project dependencies (headers & libs, CMake & CMake-gui):
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/MEPP/mepp_prebuilt_binaries_vs2010_x64_v01.rar (646 Mo)
and decompress the archive in the directory of your choice (example: C:\dev - warning, do not use an 'exotic' path like "Documents and Settings\$USER$\My ...\My ...")

2b) download MEPP project dependencies update (ffmpeg video support, assimp and zlib):
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2010/MEPP/mepp_prebuilt_binaries_vs2010_x64_addon_v02.rar (13 Mo)
and decompress the archive in the same directory as above (if necessary, answer yes to overwrite files)

3a) set 3 environment variables (menu « My computer » -> right-click -> Properties -> tab « Advanced » -> button « Environment variables » (at the bottom) -> then in « System variables » (at the bottom):
 - button « new » : add the QTDIR variable with the value:
C:\dev\qt-4.7.1 if you unzipped the above file in C:\dev
 - button « new » : add the CGAL_DIR variable with the value:
C:\dev\CGAL-3.7 if you unzipped the above file in C:\dev
 - button « change » : add in Path variable (at the end for example):
;C:\dev\qt-4.7.1\bin (pay attention to ;)

Note: you may have another program (example: MiKTeX) using a different version of Qt (and thus incompatible dlls) which will cause an error during Mepp launch.
In this case, you must change the order of your PATH variable and set C:\dev\qt-x.x.x\bin before the software in question.

3b) install Graphviz : http://download.gforge.liris.cnrs.fr/meppbin/windows/utils/ (graphviz-x.yy.z.msi)

Note: Graphviz is used by Doxygen to generate images of dependency graphs (i.e. 'dot.exe' tool).
But, if you previously installed MATLAB, it also uses a tool 'dot.exe' which will be a problem and will result in having 'empty' images.
In this case, it is necessary to change the order of your PATH variable and set Graphviz before MATLAB.

3c) reboot to taking into account of the environment variables above

4) use CMake-gui (in C:\dev\_cmake-2.8.3.20110118_\bin)
 - set the field "Where is the source code" with C:\MEPP.git (ou/or MEPP.svn)\MEPP (warning, not C:\MEPP.git (ou/or MEPP.svn)\MEPP\src !)
 - set the field "Where to build the binaries" with C:\MEPP.git (ou/or MEPP.svn)\MEPP\build
 - click on Configure (bottom left) and choose for 'Generator': "Visual Studio 10 Win64" (warning, not to be mistaken with "Visual Studio 10")
 - activate/deactivate the components you want or not (see the upper lines of the form "BUILD_component_componentname", for example: BUILD_component_Curvature)
 - click on Configure (bottom left) again
 - click on Generate (bottom left)
 - open with Visual Studio the mepp.sln generated solution which is in C:\MEPP.git (ou/or MEPP.svn)\MEPP\build then compile Mepp
 - on the "sub-project" mepp, do a "right-click" and choose "Set as startup project"
 
Note: warning, by default the project compiles in Debug mode but you can switch to Release mode if you wish.

5) Mepp documentation and your component documentation in Doxygen format are generated also via Visual Studio



* When you first start Visual Studio 2010, it is recommended to enable this option:
Menu "Tools -> Settings -> Advanced settings"

* With Visual Studio 2010 Ultimate, to compile Mepp in debug mode, you may need to disable this option:
Project « mepp » -> right-click -> Properties -> Configuration properties -> Manifest Tool -> Verbose Output -> No
