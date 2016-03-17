# MSVC - KIT VS 2015 x64 - V01

if( DEFINED ENV{MSVC_KIT_ROOT} )
	set( MSVC_KIT_ROOT $ENV{MSVC_KIT_ROOT} )
endif()
if( NOT DEFINED MSVC_KIT_ROOT )
	message(FATAL_ERROR "MSVC_KIT_ROOT not set.  Please set MSVC_KIT_ROOT.")
endif()

###

set(CGAL_DIR			${MSVC_KIT_ROOT}/CGAL-4.7)

#set(WITH_QT5			TRUE)

if(WITH_QT5)
	# with qt5
	set(QT5_DIR			${MSVC_KIT_ROOT}/Qt/Qt5.6.0/5.6/msvc2015_64)
	set(QGLVIEWERROOT	${MSVC_KIT_ROOT}/libQGLViewer/libQGLViewer-2.6.3-qt-5.6.0)
else(WITH_QT5)
	# with qt4
	set(QTDIR			${MSVC_KIT_ROOT}/Qt/qt-4.8.7-x64-msvc2015)
	set(QGLVIEWERROOT	${MSVC_KIT_ROOT}/libQGLViewer/libQGLViewer-2.6.3-qt-4.8.7)
endif(WITH_QT5)

set(XERCESC_ROOT_DIR	${MSVC_KIT_ROOT}/xerces-c-3.1.2)

### optional :

set(ASSIMP_ROOT_DIR		${MSVC_KIT_ROOT}/assimp-3.2)

set(USE_FFMPEG			TRUE)
set(FFMPEG_DIR			${MSVC_KIT_ROOT}/contrib/ffmpeg)

### for some private components :

set(ZLIB_ROOT			${MSVC_KIT_ROOT}/zlib-128)

set(JPEG_DIR			${MSVC_KIT_ROOT}/jpeg-9b)

set(IQA_DIR				${MSVC_KIT_ROOT}/iqa_1.1.2)

###