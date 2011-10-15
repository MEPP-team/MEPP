# This is the CMake script for compiling MEPP
# Martial TOLA - Juin 2011
IF(APPLE)
	cmake_minimum_required(VERSION 2.8)
else(APPLE)
    cmake_minimum_required(VERSION 2.6)
endif(APPLE)

if (COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
  IF(POLICY CMP0017)
	cmake_policy(SET CMP0017 OLD) # or NEW
  ENDIF()
endif (COMMAND cmake_policy)

# Configuration du projet
set(PRJ_NAME mepp)
project( ${PRJ_NAME} )
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/${PRJ_NAME}_deploy")

set(DOC_DIR ${CMAKE_CURRENT_BINARY_DIR}/doc)

file(
	GLOB_RECURSE
	mepp_components
	src/components/*/CMakeLists.txt
)
FOREACH (MEPP_COMPONENT ${mepp_components})
	get_filename_component (MEPP_COMPONENT_DIR "${MEPP_COMPONENT}" PATH)
	get_filename_component (MEPP_SHORT_COMPONENT_NAME "${MEPP_COMPONENT_DIR}" NAME)
	set(MEPP_COMPONENT_NAME component_${MEPP_SHORT_COMPONENT_NAME})
	
	message(STATUS "FIND MEPP_COMPONENT_DIR: ${MEPP_COMPONENT_DIR}")
	
	OPTION(BUILD_${MEPP_COMPONENT_NAME} "Build ${MEPP_COMPONENT_NAME} " ON)
	if (BUILD_${MEPP_COMPONENT_NAME})
	  ADD_SUBDIRECTORY(${MEPP_COMPONENT_DIR})
	  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/components/_enrichment_.in/polyhedron_enrichment_define.h.in ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_define_${MEPP_SHORT_COMPONENT_NAME}.h)
	  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/components/_enrichment_.in/polyhedron_enrichment_include.h.in ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_include_${MEPP_SHORT_COMPONENT_NAME}.h)
	  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/components/_enrichment_.in/polyhedron_enrichment_facet.h.in ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_facet_${MEPP_SHORT_COMPONENT_NAME}.h)
	  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/components/_enrichment_.in/polyhedron_enrichment_halfedge.h.in ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_halfedge_${MEPP_SHORT_COMPONENT_NAME}.h)
	  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/components/_enrichment_.in/polyhedron_enrichment_vertex.h.in ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_vertex_${MEPP_SHORT_COMPONENT_NAME}.h)
	  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/components/_enrichment_.in/polyhedron_enrichment_polyhedron.h.in ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_polyhedron_${MEPP_SHORT_COMPONENT_NAME}.h)  
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

Execute_process(COMMAND ${BASH} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${CAT} ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_define_ ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include/polyhedron_enrichment_define h)
Execute_process(COMMAND ${BASH} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${CAT} ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_include_ ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include/polyhedron_enrichment_include h)
Execute_process(COMMAND ${BASH} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${CAT} ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_facet_ ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include/polyhedron_enrichment_facet h)
Execute_process(COMMAND ${BASH} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${CAT} ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_halfedge_ ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include/polyhedron_enrichment_halfedge h)
Execute_process(COMMAND ${BASH} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${CAT} ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_vertex_ ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include/polyhedron_enrichment_vertex h)
Execute_process(COMMAND ${BASH} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${CAT} ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_polyhedron_ ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include/polyhedron_enrichment_polyhedron h)

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include )
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/src/mepp )

# Find CGAL
find_package(CGAL COMPONENTS)
include( ${CGAL_USE_FILE} )

# Find Qt4 itself
set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )
find_package(Qt4)

# Find OpenGL
find_package(OpenGL)
if (APPLE)	#cmake . -DCMAKE_C_FLAGS='-arch x86_64 -arch i386'
	INCLUDE_DIRECTORIES(/usr/X11/include)
endif (APPLE)

# Find GLEW (for OpenGL-1.5 and OpenGL extensions)
find_package(GLEW)

# Find QGLViewer
if (QT4_FOUND)
  include(${QT_USE_FILE})
  find_package(QGLViewer)
endif (QT4_FOUND)

if (WIN32)
	set(QTDIR $ENV{QTDIR} )
	set(DEV_DIRECTORY ${QTDIR}/.. )
endif (WIN32)

# Find XercesC
if (WIN32)
	# Hint to included XercesC resources
	SET(XERCESC_ROOT ${DEV_DIRECTORY}/contrib/xerces-c)
	SET(XML_PARSER_INSTALL_FILES ${XERCESC_ROOT}/bin/xerces-c_3_1.dll)
	SET(XERCESC_INCLUDE_DIR ${XERCESC_ROOT}/include CACHE PATH "Include path of XercesC library")
	SET(XERCESC_LIBRARY ${XERCESC_ROOT}/lib/xerces-c_3.lib CACHE FILEPATH  "XercesC library" )
endif (WIN32)
find_package(XercesC)
if (XERCESC_FOUND)
	set (XML_PARSER_INCLUDE_DIR ${XERCESC_INCLUDE_DIR})
	set (XML_PARSER_LIBRARY ${XERCESC_LIBRARY})
endif(XERCESC_FOUND)

# Find FFmpeg
find_package(FFmpeg)
if (FFMPEG_FOUND)
	ADD_DEFINITIONS("-DWITH_FFMPEG ")
endif (FFMPEG_FOUND)

if (CGAL_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND AND XERCESC_FOUND)

  find_package(Boost)

  include_directories ( ${QGLVIEWER_INCLUDE_DIR} )
  
  include_directories ( ${XML_PARSER_INCLUDE_DIR} )
  
  if(GLEW_FOUND)
    include_directories ( ${GLEW_INCLUDE_PATH} )
  else(GLEW_FOUND)
    message(STATUS "NOTICE: GLEW library is not found. VBO rendering will not be available.")
  endif(GLEW_FOUND)
  
	if (FFMPEG_FOUND)
		include_directories ( ${FFMPEG_INCLUDE_DIRS} )
	endif (FFMPEG_FOUND)

	file(
		GLOB_RECURSE
		hxx_header_files
		src/mepp/*.hxx
	)
	qt4_wrap_cpp( MOC_FILES ${hxx_header_files} )
	
	file(
		GLOB_RECURSE
		all_ui_files
		src/mepp/*.ui
	)
	qt4_wrap_ui( UI_FILES ${all_ui_files} )
	
	qt4_add_resources( RESOURCE_FILES src/mepp/mepp.qrc )
	
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
		${PRJ_NAME} MACOSX_BUNDLE WIN32
		${source_files}
		${MOC_FILES}
		${UI_FILES}
		${RESOURCE_FILES}
	)
	
	# Configuration de l'édition de liens
	# Link with Qt libraries
	target_link_libraries( ${PRJ_NAME} ${QT_LIBRARIES} )
	
	# Link with CGAL
	target_link_libraries( ${PRJ_NAME} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
	
	# Link with libQGLViewer, OpenGL
	target_link_libraries( ${PRJ_NAME} ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )
	
	# Link with XercesC
	target_link_libraries( ${PRJ_NAME} ${XML_PARSER_LIBRARY} ) 
	
	if(GLEW_FOUND)
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


IF(APPLE)
	#--------------------------------------------------------------------------------
	# Now the installation stuff : http://www.cmake.org/Wiki/BundleUtilitiesExample
	#--------------------------------------------------------------------------------
	SET(plugin_dest_dir bin)
	SET(qtconf_dest_dir bin)
	SET(APPS "\${CMAKE_INSTALL_PREFIX}/bin/${PRJ_NAME}")
	IF(APPLE)
		SET(plugin_dest_dir ${PRJ_NAME}.app/Contents/MacOS)
		SET(qtconf_dest_dir ${PRJ_NAME}.app/Contents/Resources)
		SET(APPS "\${CMAKE_INSTALL_PREFIX}/${PRJ_NAME}.app")
	ENDIF(APPLE)
	IF(WIN32)
		SET(APPS "\${CMAKE_INSTALL_PREFIX}/bin/${PRJ_NAME}.exe")
	ENDIF(WIN32)
		
	#--------------------------------------------------------------------------------
	# Install the application, on Apple, the bundle is at the root of the
	# install tree, and on other platforms it'll go into the bin directory.
	INSTALL(TARGETS ${PRJ_NAME} 
	BUNDLE DESTINATION . COMPONENT Runtime
	RUNTIME DESTINATION bin COMPONENT Runtime
	)

	#--------------------------------------------------------------------------------
	# Install needed Qt plugins by copying directories from the qt installation
	# One can cull what gets copied by using 'REGEX "..." EXCLUDE'
	#INSTALL(DIRECTORY "${QT_PLUGINS_DIR}/imageformats" DESTINATION ${plugin_dest_dir}/plugins COMPONENT Runtime)

	#--------------------------------------------------------------------------------
	# install a qt.conf file
	# this inserts some cmake code into the install script to write the file
	INSTALL(CODE "
	file(WRITE \"\${CMAKE_INSTALL_PREFIX}/${qtconf_dest_dir}/qt.conf\" \"\")
	" COMPONENT Runtime)

	# install qt_menu.nib
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
	file(GLOB_RECURSE QTPLUGINS
	\"\${CMAKE_INSTALL_PREFIX}/${plugin_dest_dir}/plugins/*${CMAKE_SHARED_LIBRARY_SUFFIX}\" \"\${CMAKE_INSTALL_PREFIX}/${plugin_dest_dir}/libcomponent_*${CMAKE_SHARED_LIBRARY_SUFFIX}\")
	include(BundleUtilities)
	fixup_bundle(\"${APPS}\" \"\${QTPLUGINS}\" \"${DIRS}\")
	" COMPONENT Runtime)

	# To Create a package, one can run "cpack -G DragNDrop CPackConfig.cmake" on Mac OS X
	# where CPackConfig.cmake is created by including CPack
	# And then there's ways to customize this as well
	#set(CPACK_BINARY_DRAGNDROP ON)
	#include(CPack)
	#--------------------------------------------------------------------------------
	# Now the installation stuff : end
	#--------------------------------------------------------------------------------
ENDIF(APPLE)

	
else (CGAL_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND AND XERCESC_FOUND)

  set(MEPP_MISSING_DEPS "")

  if (NOT CGAL_FOUND)
    set(MEPP_MISSING_DEPS "the CGAL library, ${MEPP_MISSING_DEPS}")
  endif ()

  if (NOT QT4_FOUND)
    set(MEPP_MISSING_DEPS "Qt4, ${MEPP_MISSING_DEPS}")
  endif ()

  if (NOT OPENGL_FOUND)
    set(MEPP_MISSING_DEPS "OpenGL, ${MEPP_MISSING_DEPS}")
  endif ()

  if (NOT QGLVIEWER_FOUND)
    set(MEPP_MISSING_DEPS "QGLViewer, ${MEPP_MISSING_DEPS}")
  endif ()
  
  if (NOT XERCESC_FOUND)
    set(MEPP_MISSING_DEPS "Xerces-C++, ${MEPP_MISSING_DEPS}")
  endif ()

  message(STATUS "NOTICE: This software requires ${MEPP_MISSING_DEPS}and will not be compiled.")

endif (CGAL_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND AND XERCESC_FOUND)
