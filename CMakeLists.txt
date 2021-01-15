cmake_minimum_required(VERSION 3.1)
project(swarm VERSION 0.0)

# find_package(GTest)
# find_package(Threads)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(swarm src/main.cpp src/agent.cpp)

target_include_directories(swarm PRIVATE include eigen-3.3.9)

# target_include_directories(parser PRIVATE include)


# if(GTest_FOUND)
# 	add_executable(battle_test src/battle_test.cpp src/minions.cpp src/hs_board.cpp)
# 	target_include_directories(battle_test PRIVATE ${GTEST_INCLUDE_DIRS} include)
# 	target_link_libraries(battle_test ${GTEST_LIBRARIES} Threads::Threads)
# endif()
