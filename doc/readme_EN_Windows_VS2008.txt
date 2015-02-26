///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010-2011-2012-2013-2014-2015
// CNRS-Lyon, LIRIS UMR 5205
///////////////////////////////////////////////////////////////////////////

How-to for Mepp (32 bits) under Windows with Visual Studio Express 2008 SP1 :
-----------------------------------------------------------------------------


1a) download and install DAEMON Tools Lite in order to install the following 'iso' images:
http://download.gforge.liris.cnrs.fr/meppbin/windows/utils/ (DTLitexxxx-yyyy.exe)

1b) download and install Visual Studio Express 2008 SP1:
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2008/VS2008ExpressWithSP1FRAx1504731.iso (869 Mo)

2a) download MEPP project dependencies (headers & libs, CMake & CMake-gui):
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2008/MEPP/mepp_prebuilt_binaries_vs2008_x86_v02.rar (497 Mo)
and decompress the archive in the directory of your choice (example: C:\dev32 - warning, do not use an 'exotic' path like "Documents and Settings\$USER$\My ...\My ...")

2b) download MEPP project dependencies update (ffmpeg video support):
http://download.gforge.liris.cnrs.fr/meppbin/windows/vs2008/MEPP/mepp_prebuilt_binaries_vs2008_x86_ffmpeg_addon.rar (9 Mo)
and decompress the archive in the same directory as above (if necessary, answer yes to overwrite files)

3a) set 3 environment variables (menu « My computer » -> right-click -> Properties -> tab « Advanced » -> button « Environment variables » (at the bottom) -> then in « System variables » (at the bottom):
 - button « new » : add the QTDIR variable with the value:
C:\dev32\Qt_4.6.3_x86 if you unzipped the above file in C:\dev32
 - button « new » : add the CGAL_DIR variable with the value:
C:\dev32\CGAL-3.6.1_x86 if you unzipped the above file in C:\dev32
 - button « change » : add in Path variable (at the end for example):
;C:\dev32\Qt_4.6.3_x86\bin (pay attention to ;)

Note: you may have another program (example: MiKTeX) using a different version of Qt (and thus incompatible dlls) which will cause an error during Mepp launch.
In this case, you must change the order of your PATH variable and set C:\dev32\Qt_4.6.3_x86\bin before the software in question.

3b) install Graphviz : http://download.gforge.liris.cnrs.fr/meppbin/windows/utils/ (graphviz-x.yy.z.msi)

Note: Graphviz is used by Doxygen to generate images of dependency graphs (i.e. 'dot.exe' tool).
But, if you previously installed MATLAB, it also uses a tool 'dot.exe' which will be a problem and will result in having 'empty' images.
In this case, it is necessary to change the order of your PATH variable and set Graphviz before MATLAB.

3c) reboot to taking into account of the environment variables above

4) use CMake-gui (in C:\dev32\_cmake-2.8.3.20110118_\bin)
 - set the field "Where is the source code" with C:\MEPP.git (ou/or MEPP.svn)\MEPP (warning, not C:\MEPP.git (ou/or MEPP.svn)\MEPP\src !)
 - set the field "Where to build the binaries" with C:\MEPP.git (ou/or MEPP.svn)\MEPP\build
 - click on Configure (bottom left) and choose for 'Generator': "Visual Studio 9 2008"
 - activate/deactivate the components you want or not (see the upper lines of the form "BUILD_component_componentname", for example: BUILD_component_Curvature)
 - click on Configure (bottom left) again
 - click on Generate (bottom left)
 - open with Visual Studio the mepp.sln generated solution which is in C:\MEPP.git (ou/or MEPP.svn)\MEPP\build then compile Mepp
 - on the "sub-project" mepp, do a "right-click" and choose "Set as startup project"
 
Note: warning, by default the project compiles in Debug mode but you can switch to Release mode if you wish.

5) Mepp documentation and your component documentation in Doxygen format are generated also via Visual Studio