# This is the CMake script for compiling a MEPP component
# Martial TOLA - Juin 2010
cmake_minimum_required(VERSION 2.4.5)

if (COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif (COMMAND cmake_policy)

# Configuration du projet
find_path (MEPP_COMPONENT_DIR CMakeLists.txt .)
get_filename_component (MEPP_COMPONENT_NAME "${MEPP_COMPONENT_DIR}" NAME)
message(STATUS "FIND MEPP_COMPONENT_NAME: ${MEPP_COMPONENT_NAME}")

set(COMPONENT_NAME mepp_${MEPP_COMPONENT_NAME})
project( ${COMPONENT_NAME} )
set(LIBRARY_OUTPUT_PATH ../../../..)

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

if (CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  find_package(Boost)

  include_directories ( ${QGLVIEWER_INCLUDE_DIR} )
	
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

endif (CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)
