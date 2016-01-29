
# Set variables in Cmake.  Use the user set environment variables
# If these are not set, then they will simply be ignored.


# try to find the header
FIND_PATH(XERCESC_INCLUDE_DIR xercesc/parsers/SAXParser.hpp # or maybe better ? -> xercesc/util/XercesVersion.hpp
  ${XERCESC_ROOT_DIR}/include
  $ENV{XERCESC_ROOT_DIR}/include
  /usr/include
  /usr/local/include
)

# Find the release library
FIND_LIBRARY(XERCESC_LIBRARY_RELEASE
	NAMES xerces-c xerces-c_3 xerces-c_2
	PATHS
	  ${XERCESC_ROOT_DIR}/lib
	  $ENV{XERCESC_ROOT_DIR}/lib
	  /usr/lib
	  /usr/local/lib
	DOC "The name of the xerces-c library"
)

# Find the debug library
FIND_LIBRARY(XERCESC_LIBRARY_DEBUG
	NAMES xerces-cd xerces-c_3D xerces-c_2D
	PATHS
	  ${XERCESC_ROOT_DIR}/lib
	  $ENV{XERCESC_ROOT_DIR}/lib
	  /usr/lib
	  /usr/local/lib
	DOC "The name of the xerces-c library"
)


# See if anything was found.
if(XERCESC_LIBRARY_RELEASE)
  if(XERCESC_LIBRARY_DEBUG)
    set(XERCESC_LIBRARIES_ optimized ${XERCESC_LIBRARY_RELEASE} debug ${XERCESC_LIBRARY_DEBUG})
  else()
    set(XERCESC_LIBRARIES_ ${XERCESC_LIBRARY_RELEASE})
  endif()

  set(XERCESC_LIBRARIES ${XERCESC_LIBRARIES_} CACHE FILEPATH "The xerces-c library")
endif()

IF(XERCESC_INCLUDE_DIR AND XERCESC_LIBRARIES)
   SET(XERCESC_FOUND TRUE)
ENDIF(XERCESC_INCLUDE_DIR AND XERCESC_LIBRARIES)
 	
IF (XERCESC_FOUND)
  IF (NOT XERCESC_FIND_QUIETLY)
    MESSAGE (STATUS "Found Xerces-C: ${XERCESC_LIBRARIES}")
  ENDIF (NOT XERCESC_FIND_QUIETLY)
ELSE (XERCESC_FOUND)
  IF (XERCESC_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Could not find Xerces-C")
  ENDIF (XERCESC_FIND_REQUIRED)
ENDIF (XERCESC_FOUND)
