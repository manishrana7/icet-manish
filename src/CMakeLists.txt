# Make sure CMake was invoked correctly by the user.
CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0 FATAL_ERROR)
#IF(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
#    MESSAGE(FATAL_ERROR "You should specify the root directory as the CMake source directory, not the src/ subdirectory.")
#ENDIF()


# This library depends on pybind11 framework.
ADD_SUBDIRECTORY("../3rdparty/pybind11" pybind11)

# List of source files that need to be compiled:
SET(SourceFiles
	Structure.cpp
)

PYBIND11_ADD_MODULE(example ${SourceFiles})
SET_TARGET_PROPERTIES(example PROPERTIES OUTPUT_NAME "example")

