## Windows installation

### Prerequisites

 - Update Windows 7/8/8.1/10 (services packs and Windows Update)

 - Download and install Visual Studio Express 2015 from [LIRIS host](https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/Visual%20Studio%20Express%202015/Visual%20Studio%20Express%202015%20pour%20Windows%20Desktop.rar)
   (download size 7.4 Go, installation size ~>14 Go)

### Installing dependencies

 1. Download (mandatory) the ['core' binary kit (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP1_local_vs2015_64.7z) that delivers CMake, Doxygen, Graphviz, Boost, CGAL, Xerces, Qt et libQGLViewer for `VS2015_64` (download size 511 Mo, installation size ~4.5 Go)

 2. Extract CMake 3.4.3 from 'path_to\local_vs2015_64\_utils_\cmake-3.4.3-win32-x86.zip'

 3. Set a new user environment variable 'MSVC_KIT_ROOT' to 'path_to/local_vs2015_64' (beware of the directory separator, it must be '/' here)

 4. Add ';path_to\local_vs2015_64\_bin_' to the PATH system environment variable, just after the Windows system paths, but before any other path, in order to avoid a library version conflict; beware of the ';' path separator

### Building stage

 - Get Mepp source code using your favourite Git client (first 'fork' and then 'clone') : https://github.com/MEPP-team/MEPP

 - Run cmake-gui.exe

 - Choose 'Visual Studio 14 2015 Win64' as the compiler version (beware : vs2015 is version 14 !)

 - Where is the source code = ".../MEPP"

 - Where to build the binaries = ".../MEPP/build"

 - Click "Configure"

 - Click "Generate"

 - Open ".../MEPP/build/MEPP.sln" solution with MSVC 2015, select 'Release' mode, then generate the 'ALL_BUILD' target
 
 ### Special handling to activate MSDM2
 
 - Prior to run cmake-gui.exe, change "CGAL-4.9" into "CGAL-4.7" (line 14) in the file \MEPP\cmake\msvc\kit_vs2015-64_v01.cmake.
  
 - Check "BUILD_Component_MSDM2" in the cmake-gui panel.

