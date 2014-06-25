# This is the CMake script for compiling MEPP
# Martial TOLA - June 2014
IF(APPLE) #cmake .. -DCMAKE_C_FLAGS='-arch x86_64 -arch i386'
	cmake_minimum_required(VERSION 2.8)
else(APPLE)
	cmake_minimum_required(VERSION 2.6)
endif(APPLE)

if (WITH_QT5)
	cmake_minimum_required(VERSION 2.8.9) # because 2.8.8 has a bug with QT5_USE_MODULES 
endif(WITH_QT5)

# Find includes in corresponding build directories
#set(CMAKE_INCLUDE_CURRENT_DIR ON)

if (COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
  IF(POLICY CMP0017)
	cmake_policy(SET CMP0017 NEW)
  ENDIF()
  IF(POLICY CMP0020)
	cmake_policy(SET CMP0020 NEW)
  ENDIF()
endif (COMMAND cmake_policy)

# Project configuration
set(PRJ_NAME mepp)
set(PRJ_VERSION 0.50.3) # TODO
project( ${PRJ_NAME} )
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/${PRJ_NAME}_deploy")

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(DOC_DIR ${CMAKE_CURRENT_BINARY_DIR}/doc)

#------------------------------------------------------------------------------
# Specific compiler options
#------------------------------------------------------------------------------
# Clang (export CXX=clang++ CC=clang)
#------------------------------------------------------------------------------
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qunused-arguments")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Qunused-arguments")
  message(STATUS "Clang compiler detected")
  if ( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
#    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wdocumentation")
  endif( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
endif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
#------------------------------------------------------------------------------
# MSVC
#------------------------------------------------------------------------------
# /MP for multiple compilation units (cl.exe) ; without any value, it is set automatically depending on your number of main threads
# If you are compiling on Windows, the preprocessor symbols 'min' or 'max' are defined. Do #define NOMINMAX to prevent windows.h from defining these symbols.
if(MSVC)
	SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /D \"NOMINMAX\"")
endif()
#----------------------------------------------------------------------------

file(
	GLOB_RECURSE
	mepp_components
	src/components/*/CMakeLists.txt
	src/private_components/*/CMakeLists.txt
)
FOREACH (MEPP_COMPONENT ${mepp_components})
	get_filename_component (MEPP_COMPONENT_DIR "${MEPP_COMPONENT}" PATH)
	get_filename_component (MEPP_SHORT_COMPONENT_NAME "${MEPP_COMPONENT_DIR}" NAME)
	set(MEPP_COMPONENT_NAME component_${MEPP_SHORT_COMPONENT_NAME})
	
	message(STATUS "FIND MEPP_COMPONENT_DIR: ${MEPP_COMPONENT_DIR}")
	
	OPTION(BUILD_${MEPP_COMPONENT_NAME} "Build ${MEPP_COMPONENT_NAME} " ON)
	if (BUILD_${MEPP_COMPONENT_NAME})
		ADD_SUBDIRECTORY(${MEPP_COMPONENT_DIR})
		CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/mepp/_enrichment_.in/polyhedron_enrichment_define.h.in ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_define_${MEPP_SHORT_COMPONENT_NAME}.h)
		CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/mepp/_enrichment_.in/polyhedron_enrichment_include.h.in ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_include_${MEPP_SHORT_COMPONENT_NAME}.h)
		CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/mepp/_enrichment_.in/polyhedron_enrichment_facet.h.in ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_facet_${MEPP_SHORT_COMPONENT_NAME}.h)
		CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/mepp/_enrichment_.in/polyhedron_enrichment_halfedge.h.in ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_halfedge_${MEPP_SHORT_COMPONENT_NAME}.h)
		CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/mepp/_enrichment_.in/polyhedron_enrichment_vertex.h.in ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_vertex_${MEPP_SHORT_COMPONENT_NAME}.h)
		CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/mepp/_enrichment_.in/polyhedron_enrichment_polyhedron.h.in ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_polyhedron_${MEPP_SHORT_COMPONENT_NAME}.h)  
	else (BUILD_${MEPP_COMPONENT_NAME})
		FILE( REMOVE ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_define_${MEPP_SHORT_COMPONENT_NAME}.h)
		FILE( REMOVE ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_include_${MEPP_SHORT_COMPONENT_NAME}.h)
		FILE( REMOVE ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_facet_${MEPP_SHORT_COMPONENT_NAME}.h)
		FILE( REMOVE ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_halfedge_${MEPP_SHORT_COMPONENT_NAME}.h)
		FILE( REMOVE ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_vertex_${MEPP_SHORT_COMPONENT_NAME}.h)
		FILE( REMOVE ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_polyhedron_${MEPP_SHORT_COMPONENT_NAME}.h)
	endif (BUILD_${MEPP_COMPONENT_NAME})
ENDFOREACH(MEPP_COMPONENT)

if (WIN32)
	set(CAT gnu_cat.bat)
	set(BASH "")
else (WIN32)
	set(CAT cat.sh)
	set(BASH bash)
endif (WIN32)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/mepp/mepp_config.h.in ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include/mepp_config.h)

Execute_process(COMMAND ${BASH} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${CAT} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/ ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_define_ ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include/polyhedron_enrichment_define h)
Execute_process(COMMAND ${BASH} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${CAT} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/ ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_include_ ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include/polyhedron_enrichment_include h)
Execute_process(COMMAND ${BASH} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${CAT} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/ ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_facet_ ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include/polyhedron_enrichment_facet h)
Execute_process(COMMAND ${BASH} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${CAT} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/ ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_halfedge_ ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include/polyhedron_enrichment_halfedge h)
Execute_process(COMMAND ${BASH} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${CAT} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/ ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_vertex_ ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include/polyhedron_enrichment_vertex h)
Execute_process(COMMAND ${BASH} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${CAT} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/ ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_polyhedron_ ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include/polyhedron_enrichment_polyhedron h)

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include )
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/src/mepp )

# Find CGAL libraries
FIND_PACKAGE(CGAL COMPONENTS Core)
IF(CGAL_FOUND)
	INCLUDE(${CGAL_USE_FILE})
	MESSAGE(STATUS "Found CGAL: ${CGAL_LIBRARIES}")
	MESSAGE(STATUS "Found CGAL 3RD PARTY LIBRARIES: ${CGAL_3RD_PARTY_LIBRARIES}")
ELSE(CGAL_FOUND)
	MESSAGE(FATAL_ERROR "CGAL not found. Please set CGAL_DIR.")
ENDIF(CGAL_FOUND)

# Find OpenMesh library
#if (OPENMESH_DIR)
#  set( ENV{OPENMESH_DIR} ${OPENMESH_DIR} )
#endif()
#find_package(OpenMesh)
#if(OPENMESH_FOUND)
#  include_directories(${OPENMESH_INCLUDE_DIR})
#else(OPENMESH_FOUND)
#  message(FATAL_ERROR "OpenMesh not found. Please set OPENMESH_DIR.")
#endif(OPENMESH_FOUND)

# Find QT4/QT5 libraries
if (WITH_QT5)
  # QT5 Handling
  if(DEFINED ENV{QT5_DIR})
    set( QT5_DIR $ENV{QT5_DIR} )
  endif()
  SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${QT5_DIR})

  # Instruct CMake to run moc automatically when needed
  set(CMAKE_AUTOMOC ON)

  FIND_PACKAGE(Qt5Core)
  if(${Qt5Core_FOUND})
    message(STATUS "--> Qt5Core found.")
    FIND_PACKAGE(Qt5Widgets)
    if(${Qt5Widgets_FOUND})
      message(STATUS "--> Qt5Widgets found.")
      FIND_PACKAGE(Qt5Xml)
      if(${Qt5Xml_FOUND})
	message(STATUS "--> Qt5Xml found.")
	FIND_PACKAGE(Qt5OpenGL)
	if(${Qt5OpenGL_FOUND})
	  message(STATUS "--> Qt5OpenGL found.")
	  if (WITH_QT3D)
	    FIND_PACKAGE(Qt53D)
	    if(${Qt53D_FOUND})
	      message(STATUS "--> Qt53D found.")
	    else(${Qt53D_FOUND})
	      message(FATAL_ERROR "Qt53D not found.")
	    endif(${Qt53D_FOUND})
	  endif(WITH_QT3D)
	else(${Qt5OpenGL_FOUND})
	  message(FATAL_ERROR "Qt5OpenGL not found.")
	endif(${Qt5OpenGL_FOUND})
      else(${Qt5Xml_FOUND})
	message(FATAL_ERROR "Qt5Xml not found.")
      endif(${Qt5Xml_FOUND})
    else(${Qt5Widgets_FOUND})
      message(FATAL_ERROR "Qt5Widgets not found.")
    endif(${Qt5Widgets_FOUND})
  else(${Qt5Core_FOUND})
    message(FATAL_ERROR "Qt5Core not found. Please set QT5_DIR.")
  endif(${Qt5Core_FOUND})
else(WITH_QT5)
  # QT4 Handling
  if (QTDIR)
    set( ENV{QTDIR} ${QTDIR} )
  endif()
  set( QT_USE_QTMAIN   TRUE )
  set( QT_USE_QTXML    TRUE )
  set( QT_USE_QTOPENGL  TRUE )
  FIND_PACKAGE(Qt4)
  if(QT4_FOUND)
    INCLUDE(${QT_USE_FILE})
  else(QT4_FOUND)
    message(FATAL_ERROR "Qt4 not found. Please set QTDIR.")
  endif(QT4_FOUND)
endif(WITH_QT5)

# Find OpenGL
find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})
if (APPLE)
	INCLUDE_DIRECTORIES(/usr/X11/include)
endif (APPLE) 

# Find GLEW (for OpenGL-1.5 and OpenGL extensions)
#find_package(GLEW)

# Find QGLViewer
if (QGLVIEWERROOT)
  set( ENV{QGLVIEWERROOT} ${QGLVIEWERROOT} )
endif()
find_package(QGLViewer)
if(QGLVIEWER_FOUND)
  include_directories(${QGLVIEWER_INCLUDE_DIR})
else(QGLVIEWER_FOUND)
  message(FATAL_ERROR "QGLViewer not found. Please set QGLVIEWERROOT.")
endif(QGLVIEWER_FOUND)

# Find XercesC
if (WIN32)
	if (WITH_QT5)
	else(WITH_QT5)
			SET(XERCESC_ROOT_DIR ${CGAL_DIR}/../contrib/xerces-c) #very bad: temporary, for compatibility with previous Mepp binary kits
	endif(WITH_QT5)
endif (WIN32)
find_package(XercesC)
if (XERCESC_FOUND)
	set (XML_PARSER_INCLUDE_DIR ${XERCESC_INCLUDE_DIR})
	set (XML_PARSER_LIBRARY ${XERCESC_LIBRARIES})
else(XERCESC_FOUND)
  message(FATAL_ERROR "XercesC not found. Please set XERCESC_ROOT_DIR.")
endif(XERCESC_FOUND)

# Find FFmpeg
if (FFMPEG_DIR)
  set( ENV{FFMPEG_DIR} ${FFMPEG_DIR} )
endif()
find_package(FFmpeg)
if (FFMPEG_FOUND)
  set(my_COMPILE_DEFINITIONS ${my_COMPILE_DEFINITIONS} ";WITH_FFMPEG")
else(FFMPEG_FOUND)
  message(STATUS "FFmpeg (optional) not found. Please set FFMPEG_DIR.")
endif(FFMPEG_FOUND)

# Find ASSIMP
find_package(ASSIMP)
if(ASSIMP_FOUND)
  include_directories(${ASSIMP_INCLUDE_PATH})
  set(my_COMPILE_DEFINITIONS ${my_COMPILE_DEFINITIONS} ";WITH_ASSIMP")
else(ASSIMP_FOUND)
  message(STATUS "ASSIMP (optional) not found. Please set ASSIMP_ROOT_DIR.")
endif(ASSIMP_FOUND)
  
include_directories ( ${XML_PARSER_INCLUDE_DIR} )
  
if (GLEW_FOUND)
	include_directories ( ${GLEW_INCLUDE_PATH} )
else(GLEW_FOUND)
	#    message(STATUS "NOTICE: GLEW library is not found. VBO rendering will not be available.")
endif(GLEW_FOUND)
  
if (FFMPEG_FOUND)
	include_directories ( ${FFMPEG_INCLUDE_DIRS} )
endif (FFMPEG_FOUND)

file(
	GLOB_RECURSE
	hxx_header_files
	src/mepp/*.hxx
)
	
file(
	GLOB_RECURSE
	all_ui_files
	src/mepp/*.ui
)

if (WITH_QT5)
	set( MOC_FILES "" ) # because CMAKE_AUTOMOC is ON
	qt5_wrap_ui( UI_FILES ${all_ui_files} )
	qt5_add_resources( RESOURCE_FILES src/mepp/mepp.qrc )
else(WITH_QT5) 
	qt4_wrap_cpp( MOC_FILES ${hxx_header_files} OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED )
	qt4_wrap_ui( UI_FILES ${all_ui_files} )
	qt4_add_resources( RESOURCE_FILES src/mepp/mepp.qrc )
endif(WITH_QT5)
	
# Configuration de l'exécutable
file(
	GLOB_RECURSE
	source_files
	src/mepp/*
)
	
# For Windows set the ico file with the rc file
IF(WIN32)
	SET(source_files ${source_files} ${CMAKE_CURRENT_SOURCE_DIR}/src/mepp/Pictures/${PRJ_NAME}.rc)
ENDIF(WIN32)
	
# For Apple set the icns file containing icons
IF(APPLE)
	# set how it shows up in the Info.plist file
	SET(MACOSX_BUNDLE_ICON_FILE ${PRJ_NAME}.icns) 
	# set where in the bundle to put the icns file
	SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/src/mepp/Pictures/${PRJ_NAME}.icns PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
	# include the icns file in the target
	SET(source_files ${source_files} ${CMAKE_CURRENT_SOURCE_DIR}/src/mepp/Pictures/${PRJ_NAME}.icns)
ENDIF(APPLE)
	
add_executable(
	${PRJ_NAME} MACOSX_BUNDLE
	${source_files}
	${MOC_FILES}
	${UI_FILES}
	${RESOURCE_FILES}
)

if (WITH_QT5)
  set(my_COMPILE_DEFINITIONS ${my_COMPILE_DEFINITIONS} ";WITH_QT5")
endif(WITH_QT5)

SET_TARGET_PROPERTIES(${PRJ_NAME} PROPERTIES COMPILE_DEFINITIONS "${my_COMPILE_DEFINITIONS}")
	
# Configuration de l'édition de liens
# Link with Qt libraries
if (WITH_QT5)
	QT5_USE_MODULES( ${PRJ_NAME} Widgets Xml OpenGL )
else(WITH_QT5)
	target_link_libraries( ${PRJ_NAME} ${QT_LIBRARIES} )
endif(WITH_QT5)
	
# Link with CGAL
target_link_libraries( ${PRJ_NAME} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
	
# Link with libQGLViewer, OpenGL
target_link_libraries( ${PRJ_NAME} ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

# Link with ASSIMP
if(ASSIMP_FOUND)
	target_link_libraries( ${PRJ_NAME} ${ASSIMP_LIBRARY} ) 
endif(ASSIMP_FOUND)
	
# Link with XercesC
target_link_libraries( ${PRJ_NAME} ${XML_PARSER_LIBRARY} ) 
	
# Link with OpenMesh
if (OPENMESH_FOUND)
	target_link_libraries( ${PRJ_NAME} ${OPENMESH_LIBRARIES} )
endif (OPENMESH_FOUND)
	
if (GLEW_FOUND)
	target_link_libraries( ${PRJ_NAME} ${GLEW_LIBRARIES} )
endif(GLEW_FOUND)
	
if (FFMPEG_FOUND)
	target_link_libraries( ${PRJ_NAME} ${FFMPEG_LIBRARIES} )
endif (FFMPEG_FOUND)

#Documentation
#option(BUILD_doc_${PRJ_NAME} "Build Doxygen documentation for ${PRJ_NAME}" ON)
#if(BUILD_doc_${PRJ_NAME})
	include(TargetDoc)
#endif(BUILD_doc_${PRJ_NAME})


if (WITH_QT5)
	if( NOT DEFINED QT_PLUGINS_DIR )
		set( QT_PLUGINS_DIR "${QT5_DIR}/plugins" )
	endif( NOT DEFINED QT_PLUGINS_DIR )
endif(WITH_QT5)

IF(WIN32)
 #--------------------------------------------------------------------------------
 # Copy Qt plugins to 'Debug and Release' directories and configure qt.conf file
 #--------------------------------------------------------------------------------
 file (GLOB qtplugin_dirs RELATIVE  "${QT_PLUGINS_DIR}"  "${QT_PLUGINS_DIR}/imageformats*"  "${QT_PLUGINS_DIR}/platforms*"  )
 file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/Debug/QtPlugins" )
 file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/Release/QtPlugins" )
 foreach (qtplugin ${qtplugin_dirs})
	file(COPY "${QT_PLUGINS_DIR}/${qtplugin}" DESTINATION "${CMAKE_BINARY_DIR}/Debug/QtPlugins")
	file(COPY "${QT_PLUGINS_DIR}/${qtplugin}" DESTINATION "${CMAKE_BINARY_DIR}/Release/QtPlugins")
 endforeach ()
 file(WRITE "@CMAKE_BINARY_DIR@/Debug/qt.conf" "[Paths]\nPlugins = QtPlugins")
 file(WRITE "@CMAKE_BINARY_DIR@/Release/qt.conf" "[Paths]\nPlugins = QtPlugins")
ENDIF(WIN32)


	#--------------------------------------------------------------------------------
	# Now the installation stuff : http://www.cmake.org/Wiki/BundleUtilitiesExample
	#--------------------------------------------------------------------------------
	SET(plugin_dest_dir .)
	SET(qtconf_dest_dir .)
	SET(APPS "\${CMAKE_INSTALL_PREFIX}/${PRJ_NAME}")
	IF(APPLE)
		SET(plugin_dest_dir ${PRJ_NAME}.app/Contents/MacOS)
		SET(qtconf_dest_dir ${PRJ_NAME}.app/Contents/Resources)
		SET(APPS "\${CMAKE_INSTALL_PREFIX}/${PRJ_NAME}.app")
	ENDIF(APPLE)
	IF(WIN32)
		SET(APPS "\${CMAKE_INSTALL_PREFIX}/${PRJ_NAME}.exe")
	ENDIF(WIN32)
		
	#--------------------------------------------------------------------------------
	# Install the application, on Apple, the bundle is at the root of the
	# install tree, and on other platforms it'll go to the same directory.
	INSTALL(TARGETS ${PRJ_NAME} 
		BUNDLE DESTINATION . COMPONENT ${PRJ_NAME}
		RUNTIME DESTINATION . COMPONENT ${PRJ_NAME}
	)

	#--------------------------------------------------------------------------------
	# Install needed Qt plugins by copying directories from the qt installation
	# One can cull what gets copied by using 'REGEX "..." EXCLUDE'
	INSTALL(DIRECTORY "${QT_PLUGINS_DIR}/imageformats" DESTINATION ${plugin_dest_dir}/plugins COMPONENT ${PRJ_NAME})
	if (WITH_QT5)
		INSTALL(DIRECTORY "${QT_PLUGINS_DIR}/platforms" DESTINATION ${plugin_dest_dir}/plugins COMPONENT ${PRJ_NAME})
	endif(WITH_QT5)

	# OSG plugins
	if(OPENSCENEGRAPH_FOUND)
		INCLUDE(BundleOSGPlugins)
		install_osg_plugins( osg_plugins )
		MESSAGE( STATUS "OSG plugins: " ${osg_plugins} )
	endif(OPENSCENEGRAPH_FOUND)

	#--------------------------------------------------------------------------------
	# install a qt.conf file
	# this inserts some cmake code into the install script to write the file
	INSTALL(CODE "
		file(WRITE \"\${CMAKE_INSTALL_PREFIX}/${qtconf_dest_dir}/qt.conf\" \"[Paths]\nPlugins = plugins\")
		" COMPONENT ${PRJ_NAME})
	IF(APPLE)
		INSTALL(CODE "
			file(WRITE \"\${CMAKE_INSTALL_PREFIX}/${qtconf_dest_dir}/qt.conf\" \"[Paths]\nPlugins = MacOS/plugins\")
			" COMPONENT ${PRJ_NAME})
	ENDIF(APPLE)

	# install qt_menu.nib (APPLE)
	IF(APPLE)
		IF(EXISTS ${QT_LIBRARY_DIR}/Resources/qt_menu.nib)
			INSTALL(CODE "execute_process(COMMAND cp -R \"${QT_LIBRARY_DIR}/Resources/qt_menu.nib\" \"${qtconf_dest_dir}\"
				 WORKING_DIRECTORY \$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX})")
		ELSE(EXISTS ${QT_LIBRARY_DIR}/Resources/qt_menu.nib)
			IF(EXISTS ${QT_LIBRARY_DIR}/QtGui.framework/Resources/qt_menu.nib)
				INSTALL(CODE "execute_process(COMMAND cp -R \"${QT_LIBRARY_DIR}/QtGui.framework/Resources/qt_menu.nib\" \"${qtconf_dest_dir}\"
					 WORKING_DIRECTORY \$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX})")
			ENDIF(EXISTS ${QT_LIBRARY_DIR}/QtGui.framework/Resources/qt_menu.nib)
		ENDIF(EXISTS ${QT_LIBRARY_DIR}/Resources/qt_menu.nib)
	ENDIF(APPLE)
		
	#--------------------------------------------------------------------------------
	# Use BundleUtilities to get all other dependencies for the application to work.
	# It takes a bundle or executable along with possible plugins and inspects it
	# for dependencies.  If they are not system dependencies, they are copied.

	# directories to look for dependencies
	SET(DIRS ${QT_LIBRARY_DIRS})

	# Now the work of copying dependencies into the bundle/package
	# The quotes are escaped and variables to use at install time have their $ escaped
	# An alternative is the do a configure_file() on a script and use install(SCRIPT  ...).
	# Note that the image plugins depend on QtSvg and QtXml, and it got those copied
	# over.
	INSTALL(CODE "
		file(GLOB_RECURSE PLUGINS
		\"\${CMAKE_INSTALL_PREFIX}/${plugin_dest_dir}/plugins/*${CMAKE_SHARED_LIBRARY_SUFFIX}\" \"\${CMAKE_INSTALL_PREFIX}/${plugin_dest_dir}/osgPlugins-${OPENSCENEGRAPH_VERSION}/osgdb_*${CMAKE_SHARED_LIBRARY_SUFFIX}\" \"\${CMAKE_INSTALL_PREFIX}/${plugin_dest_dir}/*component_*${CMAKE_SHARED_LIBRARY_SUFFIX}\")
		include(BundleUtilities)
		fixup_bundle(\"${APPS}\" \"\${PLUGINS}\" \"${DIRS}\")
		" COMPONENT ${PRJ_NAME})

	# CPack Configuration
	set(CPACK_PACKAGE_NAME ${PRJ_NAME})
	set(CPACK_PACKAGE_VERSION ${PRJ_VERSION})
	set(CPACK_PACKAGE_VENDOR "Liris / CNRS")
	#set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/README)
	#set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/licence_gpl_v3.txt)

	# To Create a package, one can run "cpack -G DragNDrop CPackConfig.cmake" on Mac OS X
	# where CPackConfig.cmake is created by including CPack
	# And then there's ways to customize this as well
	IF(APPLE)
		set(CPACK_BINARY_DRAGNDROP ON)
		set(CPACK_BINARY_TGZ "OFF")
		set(CPACK_BINARY_STGZ "OFF")
		set(CPACK_BINARY_PACKAGEMAKER "OFF")
		include(CPack)
	ENDIF(APPLE)

	IF(WIN32)
		set(CPACK_GENERATOR "NSIS64")
		set(CPACK_NSIS_MENU_LINKS "${PRJ_NAME}" "${PRJ_NAME}")
		#MESSAGE( STATUS "PRJ_NAME: " ${PRJ_NAME} )
		include(CPack)
	ENDIF(WIN32)
	#--------------------------------------------------------------------------------
	# Now the installation stuff : end
	#--------------------------------------------------------------------------------
