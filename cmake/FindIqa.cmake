# Locate iqa
# This module defines
# IQA_INCLUDE_DIRS, where to find the headers
# IQA_LIBRARIES iqa libraries
#
# author Martial Tola / LIRIS - 2013
# http://liris.cnrs.fr/martial.tola/


FIND_PATH(IQA_INCLUDE_DIRS iqa.h
    ${IQA_DIR}/include
    ${IQA_DIR}
    ~/Library/Frameworks
    /Library/Frameworks
    /usr/local/include
    /usr/include
    /usr/include/iqa
    /usr/local/include/iqa #brew, manual
    /sw/include # Fink
    /opt/local/include # DarwinPorts
    /opt/local/include/iqa
    /opt/csw/include # Blastwave
    /opt/include/iqa
    /usr/X11R6/include/iqa
)


find_library(IQA_LIBRARIES
             NAMES iqa
             PATHS
             /usr/lib
             /usr/local/lib
             /sw/lib
             /opt/local/lib
			 ${IQA_DIR}/lib/release
			 ${IQA_DIR}/lib
             ${IQA_DIR}			 
             ENV LD_LIBRARY_PATH
             ENV LIBRARY_PATH
             PATH_SUFFIXES iqa
            )


SET(IQA_FOUND 0)
IF(IQA_LIBRARIES AND IQA_INCLUDE_DIRS)
  SET(IQA_FOUND 1)
ENDIF(IQA_LIBRARIES AND IQA_INCLUDE_DIRS)
