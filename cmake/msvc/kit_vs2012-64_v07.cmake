# MSVC - KIT VS 2012 x64 - V07

if( DEFINED ENV{MSVC_KIT_ROOT} )
	set( MSVC_KIT_ROOT $ENV{MSVC_KIT_ROOT} )
endif()
if( NOT DEFINED MSVC_KIT_ROOT )
	message(FATAL_ERROR "MSVC_KIT_ROOT not set.  Please set MSVC_KIT_ROOT.")
endif()

###

#set(CGAL_DIR			${MSVC_KIT_ROOT}/CGAL-4.2)		# for compatibility with OLD kits... (see beginning of CGALConfig.cmake)
set(CGAL_DIR			${MSVC_KIT_ROOT}/CGAL-4.5.1)	# for compatibility with OLD kits... (see beginning of CGALConfig.cmake)

### optional :

set(USE_FFMPEG			TRUE)

###