# ---- CMake requirements ----

cmake_minimum_required(VERSION 3.11)

set(TOP_PROJECT_NAME "OmniSketch")
project(${TOP_PROJECT_NAME})

# ---- CXX standard ----

set(CMAKE_CXX_STANDARD 17)
# To disable warnings in common/toml.h on a few c++20 compiler flags
set(CMAKE_CXX_FLAGS "-Wno-unknown-warning-option -Wno-deprecated-declarations")

# ---- Python Components ----

find_package(Python COMPONENTS Interpreter)

# ---- Include testing ----

include(CTest)

# ---- Include path of header files ----

include_directories(src)

# ---- Submodules ----

add_subdirectory(third_party/eigen)
add_subdirectory(third_party/fmt)
add_subdirectory(third_party/toml)
# Configuring toml & eigen needs a bit more work
include_directories(third_party/toml/include)
include_directories(third_party/eigen)

## ============ CBCcoin ==============
include_directories(third_party/CBC/include)
include_directories(third_party/CBC/include/coin)
# ---- Boost libraries ----

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# ---- Compile static libraries ----

find_library(LIBOSICLP NAMES libOsiClp.so HINTS third_party/CBC/lib)
find_library(LIBCLP NAMES libClp.so HINTS third_party/CBC/lib)
find_library(LIBCBC NAMES libCbc.so HINTS third_party/CBC/lib)
find_library(LIBTHD NAMES libpthread.so)
add_library(OmniTools src/impl/utils.cpp src/impl/logger.cpp src/impl/data.cpp src/impl/test.cpp src/impl/hash.cpp)
target_link_libraries(OmniTools fmt ${LIBOSICLP} ${LIBCLP} ${LIBCBC} ${LIBTHD})

### add_library(OmniTools src/impl/utils.cpp src/impl/logger.cpp src/impl/data.cpp src/impl/test.cpp src/impl/hash.cpp)
### target_link_libraries(OmniTools fmt)

# ---- Add testing ----

# add_subdirectory(test)

# ---- Pcap Library ----

function(check_library)
  find_library(${ARGV0} ${ARGV1})
  if(NOT ${ARGV0})
    message(FATAL_ERROR "lib ${ARGV1} not found")
  else()
    message(STATUS "lib ${ARGV1} found at ${${ARGV0}}")
  endif()
endfunction(check_library)

if(BUILD_PCAP_PARSER)
  check_library(PCPP Pcap++)
  check_library(PACPP Packet++)
  check_library(COMPP Common++)
  check_library(PCAP pcap)
endif()

# ---- Pcap Parser ----

if(NOT PCPP_INCLUDE_PATH)
  set(PCPP_INCLUDE_PATH /usr/local/include/pcapplusplus)
  message(STATUS "pcapplusplus is included at ${PCPP_INCLUDE_PATH} by default")
else()
  message(STATUS "pcapplusplus is included at ${PCPP_INCLUDE_PATH}")
endif()

if(BUILD_PCAP_PARSER)
  add_executable(parser ${CMAKE_CURRENT_SOURCE_DIR}/src/pcap_parser/parser.cpp)
  target_include_directories(parser PRIVATE ${PCPP_INCLUDE_PATH})
  target_link_libraries(parser ${PCPP} ${PACPP} ${COMPP} ${PCAP} fmt OmniTools)
endif()

# ---- User-defined sketches ----

function(add_user_sketch)
  add_custom_target(${ARGV1}
  DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/src/driver/${ARGV1}Driver.cpp
  )
  add_custom_command(
    OUTPUT
    ${CMAKE_CURRENT_SOURCE_DIR}/src/driver/${ARGV1}Driver.cpp
    COMMAND
        ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/src/generate_driver.py sketch_test/${ARGV1}Test.h
    DEPENDS
        ${CMAKE_CURRENT_SOURCE_DIR}/src/sketch_test/${ARGV1}Test.h
    WORKING_DIRECTORY
        ${CMAKE_CURRENT_SOURCE_DIR}/src
  )
  
  add_executable(${ARGV0} ${CMAKE_CURRENT_SOURCE_DIR}/src/driver/${ARGV1}Driver.cpp)
  target_link_libraries(${ARGV0} OmniTools)
endfunction(add_user_sketch)

# Bloom Filter
add_user_sketch(BF BloomFilter)

# Count Min Sketch
add_user_sketch(CM CMSketch)

# CH-optimized Count Min Sketch
add_user_sketch(CHCM CHCMSketch)

# THD_CH-optimized Count Min Sketch
add_user_sketch(TCHCM THD_CHCMSketch)

# CU Sketch
add_user_sketch(CU CUSketch)

# CH-optimized CU Sketch
add_user_sketch(CHCU CHCUSketch)

# HashPipe
add_user_sketch(HP HashPipe)

# CH-optimized HashPipe
add_user_sketch(CHHP CHHashPipe)

# Count Sketch
add_user_sketch(CS CountSketch)

# CH-optimized Count Sketch
add_user_sketch(CHCS CHCountSketch)

# THD CH-optimized Count Sketch
add_user_sketch(TCHCS THD_CHCountSketch)

# LD Sketch
add_user_sketch(LD LDSketch)

# Flow Radar
add_user_sketch(FR FlowRadar)

# CH-optimized Flow Radar
add_user_sketch(CHFR CHFlowRadar)

# THD CH-optimized Flow Radar
add_user_sketch(TCHFR THD_CHFlowRadar)

# Counting Bloom Filter
add_user_sketch(CBF CountingBloomFilter)

# CH-optimized Counting Bloom Filter
add_user_sketch(CHCBF CHCountingBloomFilter)

# Querying Counting Bloom Filter
add_user_sketch(QCBF QueryingCountingBloomFilter)

# CH-optimized Querying Counting Bloom Filter
add_user_sketch(CHQCBF CHQueryingCountingBloomFilter)

# NitroSketch
add_user_sketch(NS NitroSketch)

# CH-optimized NitroSketch
add_user_sketch(CHNS CHNitroSketch)

# THD CH-optimized NitroSketch
add_user_sketch(TCHNS THD_CHNitroSketch)

# Deltoid
add_user_sketch(DT Deltoid)

# CH-optimized Deltoid
add_user_sketch(CHDT CHDeltoid)

# Deltoid
add_user_sketch(DT2 Deltoid2Tuple)

# CH-optimized Deltoid
add_user_sketch(CHDT2 CHDeltoid2Tuple)

# THD CH-optimized Deltoid
add_user_sketch(TCHDT THD_CHDeltoid)

# DH Sketch
add_user_sketch(DHS DHSketch)

# HeavyKeeper
add_user_sketch(HK HeavyKeeper)

# CH-optimized Heavy Keeper
add_user_sketch(CHHK CHHeavyKeeper)

# WavingSketch
add_user_sketch(WS WavingSketch)

# CH-optimized Waving Sketch
add_user_sketch(CHWS CHWavingSketch)

# CounterTree
add_user_sketch(CT CounterTree)

# UnivMon
add_user_sketch(UM UnivMon)

# CH-optimized Sketch
add_user_sketch(CHUM CHUnivMon)

# Heavy Hitter UnivMon
add_user_sketch(HHUM HHUnivMon)

# CH-optimized Sketch
add_user_sketch(CHHHUM CHHHUnivMon)

# THD CH-optimized Sketch
add_user_sketch(TCHHHUM THD_CHHHUnivMon)

# CounterBraids
add_user_sketch(CB CounterBraids)

# MV Sketch
add_user_sketch(MV MVSketch)

# CH-optimized Sketch
add_user_sketch(CHMV CHMVSketch)

# PR Sketch
add_user_sketch(PR PRSketch)

# CH-optimized PR Sketch
add_user_sketch(CHPR CHPRSketch)

# NZE Sketch
add_user_sketch(NZE NZESketch)

# CH-optimized NZE Sketch
add_user_sketch(CHNZE CHNZESketch)

# Elastic Sketch
add_user_sketch(ES ElasticSketch)

# CH-optimized Elastic Sketch
add_user_sketch(CHES CHElasticSketch)

# PCMSketch
add_user_sketch(PCM PCMSketch)

# PCUSketch
add_user_sketch(PCU PCUSketch)

# Sketch Learn
add_user_sketch(SL SketchLearn)

# CH-optimized Sketch Learn
add_user_sketch(CHSL CHSketchLearn)

# THD CH-optimized Sketch Learn
add_user_sketch(TCHSL THD_CHSketchLearn)

# Sketch Learn
add_user_sketch(SL2 SketchLearn2Tuple)

# CH-optimized Sketch Learn
add_user_sketch(CHSL2 CHSketchLearn2Tuple)

# SALSA CU
add_user_sketch(SSCU SALSACU)

# SALSA CM
add_user_sketch(SSCM SALSACM)
