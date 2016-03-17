# MSVC - KIT VS 2015 x64 - V01

if( DEFINED ENV{MSVC_KIT_ROOT} )
	set( MSVC_KIT_ROOT $ENV{MSVC_KIT_ROOT} )
endif()
if( NOT DEFINED MSVC_KIT_ROOT )
	message(FATAL_ERROR "MSVC_KIT_ROOT not set.  Please set MSVC_KIT_ROOT.")
endif()

###

set(CGAL_DIR			${MSVC_KIT_ROOT}/CGAL-4.7)

set(WITH_QT5			TRUE)
set(QGLVIEWERROOT		${MSVC_KIT_ROOT}/libQGLViewer/libQGLViewer-2.6.3-qt-5.4.2)

set(XERCESC_ROOT_DIR	${MSVC_KIT_ROOT}/xerces-c-3.1.2)

### optional :

set(ASSIMP_ROOT_DIR		${MSVC_KIT_ROOT}/assimp-3.2)

set(USE_FFMPEG			TRUE)
set(FFMPEG_DIR			${MSVC_KIT_ROOT}/contrib/ffmpeg)

### for some private components :

#set(ZLIB_ROOT			${MSVC_KIT_ROOT}/zlib-128)

#set(JPEG_DIR			${MSVC_KIT_ROOT}/jpeg-9b)

#set(IQA_DIR				${MSVC_KIT_ROOT}/iqa_1.1.2)

###