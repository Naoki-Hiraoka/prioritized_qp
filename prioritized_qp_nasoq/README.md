## Settings required in CMakeLists.txt

```
find_package(catkin REQUIRED COMPONENTS
  prioritized_qp)

find_package(OsqpEigen REQUIRED)

include_directories(
  ${catkin_INCLUDE_DIRS}
  )

target_link_libraries([your target]
  ${catkin_LIBRARIES}
  OsqpEigen::OsqpEigen
  )
```