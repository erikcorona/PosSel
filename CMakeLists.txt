cmake_minimum_required(VERSION 3.4)
project(PosSel)

find_package(Boost COMPONENTS system filesystem thread REQUIRED)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -std=c++1y -march=native")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -std=c++1y -Wall -O1")

set(SOURCE_FILES main.cpp)
add_executable(posSel ${SOURCE_FILES})

include_directories(sel)

add_subdirectory(sel)
add_subdirectory(sel_tests)

target_link_libraries(posSel
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${Boost_LIBRARIES}
        sel2
        )