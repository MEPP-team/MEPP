# This is the CMake script for compiling a MEPP component
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

# Configuration du projet
find_path (MEPP_COMPONENT_DIR CMakeLists.txt .)
get_filename_component (MEPP_COMPONENT_NAME "${MEPP_COMPONENT_DIR}" NAME)
message(STATUS "FIND MEPP_COMPONENT_NAME: ${MEPP_COMPONENT_NAME}")

set(COMPONENT_NAME component_${MEPP_COMPONENT_NAME})
project( ${COMPONENT_NAME} )

# cf. CMakeLists de MEPP
SET(plugin_dest_dir "")
IF(APPLE)
	SET(plugin_dest_dir ${PRJ_NAME}.app/Contents/MacOS)
ENDIF(APPLE)

set(LIBRARY_OUTPUT_PATH ../../../../${plugin_dest_dir})

INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include )
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

if (FFMPEG_FOUND)
	include_directories ( ${FFMPEG_INCLUDE_DIRS} )
endif (FFMPEG_FOUND)

file(
	GLOB_RECURSE
	hxx_header_files
	src/*.hxx
)

file(
	GLOB_RECURSE
	all_ui_files
	src/*.ui
)

if (WITH_QT5)
	set( MOC_FILES "" ) # because CMAKE_AUTOMOC is ON
	qt5_wrap_ui( UI_FILES ${all_ui_files} ../../../mepp/mainwindow.ui )
else(WITH_QT5) 
	qt4_wrap_cpp( MOC_FILES ${hxx_header_files} ../../../mepp/mepp_action.hxx OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED )
	qt4_wrap_ui( UI_FILES ${all_ui_files} ../../../mepp/mainwindow.ui )
endif(WITH_QT5)

# Utilisation d'autre(s) composant(s) (sans la gestion des boîtes de dialogue pour le moment)
include( ./cmake/use_components.txt )

# Configuration de la librairie
file(
	GLOB_RECURSE
	source_files
	src/*
)
set( source_files ${source_files} ../../../mepp/mepp_action.hxx )

add_library(
	${COMPONENT_NAME}
	SHARED
	${source_files}
	${use_components}
	${MOC_FILES}
	${UI_FILES}
)

if (WITH_QT5)
  set(my_COMPILE_DEFINITIONS ${my_COMPILE_DEFINITIONS} ";WITH_QT5")
endif(WITH_QT5)

SET_TARGET_PROPERTIES(${COMPONENT_NAME} PROPERTIES COMPILE_DEFINITIONS "${my_COMPILE_DEFINITIONS}")

# Configuration de l'édition de liens
# Link with Qt libraries
if (WITH_QT5)
	QT5_USE_MODULES( ${COMPONENT_NAME} Widgets Xml OpenGL )
else(WITH_QT5)
	target_link_libraries( ${COMPONENT_NAME} ${QT_LIBRARIES} )
endif(WITH_QT5)

# Link with CGAL
target_link_libraries( ${COMPONENT_NAME} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

# Link with libQGLViewer, OpenGL
target_link_libraries( ${COMPONENT_NAME} ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

# Link with ASSIMP
if(ASSIMP_FOUND)
	target_link_libraries( ${COMPONENT_NAME} ${ASSIMP_LIBRARY} ) 
endif(ASSIMP_FOUND)

# Link with XercesC
target_link_libraries( ${COMPONENT_NAME} ${XML_PARSER_LIBRARY} )

if (FFMPEG_FOUND)
	target_link_libraries( ${COMPONENT_NAME} ${FFMPEG_LIBRARIES} )
endif (FFMPEG_FOUND) 

#Documentation
#option(BUILD_doc_${COMPONENT_NAME} "Build Doxygen documentation for ${COMPONENT_NAME}" ON)
#if(BUILD_doc_${COMPONENT_NAME})
  include(TargetDoc)
#endif(BUILD_doc_${COMPONENT_NAME})

# Install the library
IF(NOT APPLE)
	INSTALL(TARGETS ${COMPONENT_NAME} 
		LIBRARY DESTINATION . COMPONENT ${COMPONENT_NAME}
		RUNTIME DESTINATION . COMPONENT ${COMPONENT_NAME}
	)

	if (WIN32) # Release
		SET(APPS "\${CMAKE_INSTALL_PREFIX}/${PRJ_NAME}.exe")
		SET(DIRS ${QT_LIBRARY_DIRS}) # directories to look for dependencies
		INSTALL(CODE "
		   include(BundleUtilities)
		   fixup_bundle(\"${APPS}\" \"${CMAKE_BINARY_DIR}/Release/${COMPONENT_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}\" \"${DIRS}\")
		   " COMPONENT ${PRJ_NAME})
	endif (WIN32)
ENDIF(NOT APPLE)
