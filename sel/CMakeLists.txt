cmake_minimum_required(VERSION 3.4)
project(PosSel2)

add_library(sel2 SequencesAlgos.hxx Sequences.hxx WindowAlgos.hxx dummy.cxx)

set_target_properties(sel2 PROPERTIES LINKER_LANGUAGE CXX)