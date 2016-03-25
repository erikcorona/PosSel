cmake_minimum_required(VERSION 3.4)

include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(sequenceAlgos_tests sequenceAlgos_tests.cxx)

target_link_libraries(sequenceAlgos_tests gtest gtest_main)
target_link_libraries(sequenceAlgos_tests sel2)