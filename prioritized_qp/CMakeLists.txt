cmake_minimum_required(VERSION 3.0.2)
project(prioritized_qp)

add_compile_options(-std=c++11)
set(CMAKE_BUILD_TYPE Release)

find_package(catkin REQUIRED COMPONENTS prioritized_qp_base)
find_package(OsqpEigen REQUIRED)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES prioritized_qp
  CATKIN_DEPENDS prioritized_qp_base
  #DEPENDS OsqpEigen
)

###########
## Build ##
###########

## Declare a C++ library
add_library(${PROJECT_NAME}
  src/PrioritizedQPSolver.cpp
  )

target_link_libraries(${PROJECT_NAME} PUBLIC
  ${catkin_LIBRARIES}
  OsqpEigen::OsqpEigen # -isystemでこのワークスペースのdevel/includeがinclude_directoriesに入るのだが、優先順位末尾に入ってしまう
  #$<TARGET_PROPERTY:OsqpEigen::OsqpEigen,INTERFACE_LINK_LIBRARIES>
  )

target_include_directories(${PROJECT_NAME} PUBLIC
  include
  ${catkin_INCLUDE_DIRS}
  #$<TARGET_PROPERTY:OsqpEigen::OsqpEigen,INTERFACE_INCLUDE_DIRECTORIES>
)

#############
## Install ##
#############

## Mark libraries for installation
## See http://docs.ros.org/melodic/api/catkin/html/howto/format1/building_libraries.html
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

#############
## Testing ##
#############

## Add gtest based cpp test target and link libraries
# catkin_add_gtest(${PROJECT_NAME}-test test/test_prioritized_qp.cpp)
# if(TARGET ${PROJECT_NAME}-test)
#   target_link_libraries(${PROJECT_NAME}-test ${PROJECT_NAME})
# endif()

## Add folders to be run by python nosetests
# catkin_add_nosetests(test)
