# MSVC - KIT VS 2010 x64 - V01

if( DEFINED ENV{MSVC_KIT_ROOT} )
	set( MSVC_KIT_ROOT $ENV{MSVC_KIT_ROOT} )
endif()
if( NOT DEFINED MSVC_KIT_ROOT )
	message(FATAL_ERROR "MSVC_KIT_ROOT not set.  Please set MSVC_KIT_ROOT.")
endif()

###

set(CGAL_DIR					${MSVC_KIT_ROOT}/CGAL-3.7)					# for compatibility with OLD kits... (see beginning of CGALConfig.cmake)

set(XERCESC_ROOT_DIR	${MSVC_KIT_ROOT}/contrib/xerces-c)	# for compatibility with OLD kits... (see beginning of CGALConfig.cmake)

### optional :

set(USE_FFMPEG				TRUE)

###