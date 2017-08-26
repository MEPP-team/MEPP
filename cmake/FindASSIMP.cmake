#
# Try to find ASSIMP library and include path.
# Once done this will define
#
# ASSIMP_FOUND
# ASSIMP_INCLUDE_DIR
# ASSIMP_LIBRARY
#

IF (WIN32 AND NOT "${CMAKE_SIZEOF_VOID_P}" STREQUAL "4")
	FIND_PATH( ASSIMP_INCLUDE_DIR assimp/scene.h
		${ASSIMP_ROOT_DIR}/include
		DOC "The directory where assimp/scene.h resides")

    FIND_LIBRARY( ASSIMP_LIBRARY
        NAMES assimp ASSIMP assimp
        PATHS
        ${ASSIMP_ROOT_DIR}/lib/x32
        DOC "The ASSIMP library")
ELSE (WIN32 AND NOT "${CMAKE_SIZEOF_VOID_P}" STREQUAL "4")
	FIND_PATH( ASSIMP_INCLUDE_DIR assimp/scene.h
		/usr/include
		/usr/local/include
		/sw/include
		/opt/local/include
		$ENV{ASSIMP_ROOT_DIR}/include
		DOC "The directory where assimp/scene.h resides")

	FIND_LIBRARY( ASSIMP_LIBRARY
		NAMES ASSIMP assimp
		PATHS
		/usr/lib64
		/usr/lib
		/usr/local/lib64
		/usr/local/lib
		/sw/lib
		/opt/local/lib
		$ENV{ASSIMP_ROOT_DIR}/lib/x64
		DOC "The ASSIMP library")
ENDIF (WIN32 AND NOT "${CMAKE_SIZEOF_VOID_P}" STREQUAL "4")

SET(ASSIMP_FOUND "NO")
IF (ASSIMP_INCLUDE_DIR AND ASSIMP_LIBRARY)
	SET(ASSIMP_LIBRARIES ${ASSIMP_LIBRARY})
	SET(ASSIMP_FOUND "YES")
	MESSAGE(STATUS "Found ASSIMP: ${ASSIMP_LIBRARY}")
ELSE (ASSIMP_INCLUDE_DIR AND ASSIMP_LIBRARY)
	MESSAGE(STATUS "Could not find Assimp")
ENDIF (ASSIMP_INCLUDE_DIR AND ASSIMP_LIBRARY)
