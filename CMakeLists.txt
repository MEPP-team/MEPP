# This is the CMake script for compiling MEPP
# Martial TOLA - Juin 2010
cmake_minimum_required(VERSION 2.4.5)

if (COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif (COMMAND cmake_policy)

# Configuration du projet
set(PRJ_NAME mepp)
project( ${PRJ_NAME} )

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
	endif (BUILD_${MEPP_COMPONENT_NAME})
ENDFOREACH(MEPP_COMPONENT)

if (WIN32)
	set(CAT gnu_cat.bat)
	set(BASH "")
else (WIN32)
	set(CAT cat.sh)
	set(BASH bash)
endif (WIN32)

Execute_process(COMMAND ${BASH} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${CAT} ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_define_ ${CMAKE_CURRENT_SOURCE_DIR}/src/mepp/Polyhedron/polyhedron_enrichment_define h)
Execute_process(COMMAND ${BASH} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${CAT} ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_include_ ${CMAKE_CURRENT_SOURCE_DIR}/src/mepp/Polyhedron/polyhedron_enrichment_include h)
Execute_process(COMMAND ${BASH} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${CAT} ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_facet_ ${CMAKE_CURRENT_SOURCE_DIR}/src/mepp/Polyhedron/polyhedron_enrichment_facet h)
Execute_process(COMMAND ${BASH} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${CAT} ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_halfedge_ ${CMAKE_CURRENT_SOURCE_DIR}/src/mepp/Polyhedron/polyhedron_enrichment_halfedge h)
Execute_process(COMMAND ${BASH} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${CAT} ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_vertex_ ${CMAKE_CURRENT_SOURCE_DIR}/src/mepp/Polyhedron/polyhedron_enrichment_vertex h)
Execute_process(COMMAND ${BASH} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${CAT} ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/h.tmp/polyhedron_enrichment_polyhedron_ ${CMAKE_CURRENT_SOURCE_DIR}/src/mepp/Polyhedron/polyhedron_enrichment_polyhedron h)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/mepp/mepp_config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/mepp/mepp_config.h)

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
	add_executable(
		${PRJ_NAME}
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

  message(STATUS "NOTICE: This software requires ${MEPP_MISSING_DEPS}and will not be compiled.")

endif (CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)
