cmake_minimum_required(VERSION 3.15)
project(geo_boost)

set(CMAKE_CXX_STANDARD 17)

# Find Boost and include necessary libraries
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Find SFML and include necessary libraries
find_package(SFML 2.5 COMPONENTS graphics window system REQUIRED)
include_directories(${SFML_INCLUDE_DIRS})

# Specify the source directory
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Add the executable using the source files in SRC_DIR
add_executable(geo_boost ${SRC_DIR}/main.cpp)

# Link Boost and SFML libraries
target_link_libraries(geo_boost ${Boost_LIBRARIES} sfml-graphics sfml-window sfml-system)




