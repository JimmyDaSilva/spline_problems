cmake_minimum_required(VERSION 2.8.3)
project(spline_problems)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_BUILD_TYPE Release)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
  eigen_conversions
  cmake_modules
)

find_package(Eigen REQUIRED)

# catkin_python_setup()

catkin_package(
  INCLUDE_DIRS include
#  LIBRARIES spline_problems
#  CATKIN_DEPENDS other_catkin_pkg
  DEPENDS Eigen
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${Eigen_INCLUDE_DIRS}
)

## Declare a C++ library
# add_library(${PROJECT_NAME}
#   src/${PROJECT_NAME}/spline_problems.cpp
# )

# add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

add_executable(min_jerk src/min_jerk.cpp)

# add_dependencies(${PROJECT_NAME}_node ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

## Specify libraries to link a library or executable target against
target_link_libraries(min_jerk ${catkin_LIBRARIES} -lpython2.7)