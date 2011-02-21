# This is the CMake script for compiling a MEPP component
# Martial TOLA - Janvier 2011
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
find_path (MEPP_COMPONENT_DIR CMakeLists.txt .)
get_filename_component (MEPP_COMPONENT_NAME "${MEPP_COMPONENT_DIR}" NAME)
message(STATUS "FIND MEPP_COMPONENT_NAME: ${MEPP_COMPONENT_NAME}")

set(COMPONENT_NAME component_${MEPP_COMPONENT_NAME})
project( ${COMPONENT_NAME} )

# cf. CMakeLists de mepp
SET(plugin_dest_dir "")
IF(APPLE)
	SET(plugin_dest_dir ${PRJ_NAME}.app/Contents/MacOS)
ENDIF(APPLE)

set(LIBRARY_OUTPUT_PATH ../../../../${plugin_dest_dir})

INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include )
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/src/mepp )

# Find CGAL and CGAL Qt4
find_package(CGAL COMPONENTS Qt4)
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

if (CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND AND XERCESC_FOUND)

  find_package(Boost)

  include_directories ( ${QGLVIEWER_INCLUDE_DIR} )
  
  include_directories ( ${XML_PARSER_INCLUDE_DIR} )
	
	file(
		GLOB_RECURSE
		hxx_header_files
		src/*.hxx
	)
	qt4_wrap_cpp( MOC_FILES ${hxx_header_files} ../../../mepp/mepp_action.hxx )
	
	file(
		GLOB_RECURSE
		all_ui_files
		src/*.ui
	)
	qt4_wrap_ui( UI_FILES ${all_ui_files} ../../../mepp/mainwindow.ui )
	
	# Utilisation d'autre(s) composant(s) (sans la gestion des boîtes de dialogue pour le moment)
	include( ./cmake/use_components.txt )
	
	# Configuration de l'exécutable
	file(
		GLOB_RECURSE
		source_files
		src/*
	)
	add_library(
		${COMPONENT_NAME}
		SHARED
		${source_files}
		${use_components}
		${MOC_FILES}
		${UI_FILES}
	)
	
	# Configuration de l'édition de liens
	# Link with Qt libraries
	target_link_libraries( ${COMPONENT_NAME} ${QT_LIBRARIES} )
	
	# Link with CGAL
	target_link_libraries( ${COMPONENT_NAME} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
	
	# Link with libQGLViewer, OpenGL
	target_link_libraries( ${COMPONENT_NAME} ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )
	
	# Link with XercesC
	target_link_libraries( ${COMPONENT_NAME} ${XML_PARSER_LIBRARY} ) 
	
	#Documentation
	option(BUILD_${COMPONENT_NAME}_DOC "Build Doxygen documentation for ${COMPONENT_NAME}" ON)
	if(BUILD_${COMPONENT_NAME}_DOC)
	  include(TargetDoc)
	endif(BUILD_${COMPONENT_NAME}_DOC)

else (CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  set(MEPP_MISSING_DEPS "")

  if (NOT CGAL_Qt4_FOUND)
    set(MEPP_MISSING_DEPS "the CGAL Qt4 library, ${MEPP_MISSING_DEPS}")
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

  message(STATUS "NOTICE: This software component requires ${MEPP_MISSING_DEPS}and will not be compiled.")

endif (CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND AND XERCESC_FOUND)
