cmake_minimum_required(VERSION 3.0.2)
project(prioritized_qp_qpswift)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Release)

find_package(catkin REQUIRED COMPONENTS prioritized_qp_base qpswifteigen)
find_package(qpSWIFT REQUIRED)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES prioritized_qp_qpswift
  CATKIN_DEPENDS prioritized_qp_base qpswifteigen
)

###########
## Build ##
###########

## Declare a C++ library
add_library(${PROJECT_NAME}
  src/PrioritizedQPQPSWIFTSolver.cpp
  )

target_link_libraries(${PROJECT_NAME} PUBLIC
  ${catkin_LIBRARIES}
  qpSWIFT::qpSWIFT-shared
  )

target_include_directories(${PROJECT_NAME} PUBLIC
  include
  ${catkin_INCLUDE_DIRS}
)

#############
## Install ##
#############

install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_GLOBAL_BIN_DESTINATION}
  )

## Mark cpp header files for installation
install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
  FILES_MATCHING PATTERN "*.h"
  PATTERN ".svn" EXCLUDE
  )
