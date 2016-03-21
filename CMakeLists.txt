cmake_minimum_required(VERSION 3.4)
project(posSel)

find_package(Boost COMPONENTS system filesystem thread REQUIRED)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -std=c++1y -march=native")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -std=c++1y -Wall")

set(SOURCE_FILES main.cpp Sequences.hxx)
add_executable(posSel ${SOURCE_FILES})

target_link_libraries(posSel
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${Boost_LIBRARIES}
        )