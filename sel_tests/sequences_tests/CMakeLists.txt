cmake_minimum_required(VERSION 3.4)

include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(sequences_tests sequences_tests.cxx)
target_link_libraries(sequences_tests gtest gtest_main)
target_link_libraries(sequences_tests sel2)