cmake_minimum_required(VERSION 3.0.2)
project(prioritized_qp_base)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Release)

find_package(catkin REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(EIGEN eigen3)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES prioritized_qp_base
  DEPENDS EIGEN
)

###########
## Build ##
###########

## Declare a C++ library
add_library(${PROJECT_NAME}
  src/PrioritizedQPBaseSolver.cpp
  )

target_link_libraries(${PROJECT_NAME} PUBLIC
  ${EIGEN_LIBRARIES}
  )

target_include_directories(${PROJECT_NAME} PUBLIC
  include
  ${EIGEN_INCLUDE_DIRS}
)

#############
## Install ##
#############

install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_GLOBAL_BIN_DESTINATION}
  )

install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
  FILES_MATCHING PATTERN "*.h"
  PATTERN ".svn" EXCLUDE
  )
