cmake_minimum_required(VERSION 3.0.2)
project(prioritized_qp_proxsuite)

add_compile_options(-std=c++14 -Wno-unused-result)
set(CMAKE_BUILD_TYPE Release)

find_package(proxsuite)

find_package(catkin REQUIRED COMPONENTS
  prioritized_qp_base)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES prioritized_qp_proxsuite
  CATKIN_DEPENDS prioritized_qp_base
)

add_library(${PROJECT_NAME}
  src/prioritized_qp_proxsuite.cpp
  )

target_link_libraries(${PROJECT_NAME} PUBLIC
  ${catkin_LIBRARIES}
  )
target_include_directories(${PROJECT_NAME} PUBLIC
  include
  ${catkin_INCLUDE_DIRS}
  ${proxsuite_INCLUDE_DIRS}
  )

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
