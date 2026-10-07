cmake_minimum_required(VERSION 3.24)

project(calorimeter_waveform_studies VERSION 0.1 LANGUAGES CXX)

set(headers
    include/HCal.hpp
)

find_package(ROOT REQUIRED COMPONENTS Core RIO)

add_library(waveformLib SHARED ${headers})
target_include_directories(waveformLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(waveformLib PUBLIC ROOT::Core ROOT::RIO)

root_generate_dictionary(G__waveformLib
    ${headers}
    MODULE waveformLib
    LINKDEF LinkDef.h 
)