cmake_minimum_required(VERSION 3.1)
project(swarm VERSION 0.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(swarm src/main.cpp src/agent.cpp)

target_include_directories(swarm PRIVATE include)
target_link_libraries (swarm Eigen3::Eigen)
