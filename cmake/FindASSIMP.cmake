#
# Try to find ASSIMP library and include path.
# Once done this will define
#
# ASSIMP_FOUND
# ASSIMP_INCLUDE_PATH
# ASSIMP_LIBRARY
#

IF (WIN32)
	FIND_PATH( ASSIMP_INCLUDE_PATH assimp/scene.h
		${ASSIMP_ROOT_DIR}/include
		DOC "The directory where assimp/scene.h resides")

    FIND_LIBRARY( ASSIMP_LIBRARY
        NAMES assimp ASSIMP assimp
        PATHS
        ${ASSIMP_ROOT_DIR}/lib
        DOC "The ASSIMP library")
ELSE (WIN32)
	FIND_PATH( ASSIMP_INCLUDE_PATH assimp/scene.h
		/usr/include
		/usr/local/include
		/sw/include
		/opt/local/include
		${ASSIMP_ROOT_DIR}/include
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
		${ASSIMP_ROOT_DIR}/lib
		DOC "The ASSIMP library")
ENDIF (WIN32)

SET(ASSIMP_FOUND "NO")
IF (ASSIMP_INCLUDE_PATH AND ASSIMP_LIBRARY)
	SET(ASSIMP_LIBRARIES ${ASSIMP_LIBRARY})
	SET(ASSIMP_FOUND "YES")
	MESSAGE(STATUS "Found ASSIMP: ${ASSIMP_LIBRARY}")
ENDIF (ASSIMP_INCLUDE_PATH AND ASSIMP_LIBRARY)
