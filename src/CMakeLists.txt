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
  PyBinding.cpp
  Neighborlist.cpp
  ManybodyNeighborlist.cpp
  Cluster.cpp
  PermutationMap.cpp
  ClusterCounts.cpp
  Orbit.cpp
  OrbitList.cpp
  LocalOrbitlistGenerator.cpp
  ClusterSpace.cpp

)
# Find Eigen3
FIND_PACKAGE(Eigen3 REQUIRED)
IF(NOT Eigen3_FOUND)
  MESSAGE(FATAL_ERROR "Eigen3 not found.")
ENDIF()
include_directories(${EIGEN3_INCLUDE_DIR})

#target_include_directories(_icetdev PRIVATE ${EIGEN3_INCLUDE_DIR})
PYBIND11_ADD_MODULE(_icetdev ${SourceFiles})
SET_TARGET_PROPERTIES(_icetdev PROPERTIES OUTPUT_NAME "_icetdev")


