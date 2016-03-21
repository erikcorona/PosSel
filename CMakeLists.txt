cmake_minimum_required(VERSION 3.4)
project(posSel)

find_package(Boost COMPONENTS system filesystem thread REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")

set(SOURCE_FILES main.cpp)
add_executable(posSel ${SOURCE_FILES})

target_link_libraries(posSel
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${Boost_LIBRARIES}
        )