cmake_minimum_required(VERSION 3.10)

project(CopcVerify)

set(LAZPERF_DIR ${PROJECT_SOURCE_DIR}/lazperf)
file(GLOB_RECURSE LAZPERF_SRCS "${LAZPERF_DIR}/*.cpp")

set(SRCS CopcVerify.cpp ${LAZPERF_SRCS})

add_executable(copcverify ${SRCS})
#target_compile_options(copcverify PRIVATE -fsanitize=undefined)
#target_link_options(copcverify PRIVATE -fsanitize=undefined)

target_include_directories(copcverify PRIVATE .)
set_property(TARGET copcverify PROPERTY CXX_STANDARD 17)
set_property(TARGET copcverify PROPERTY CXX_STANDARD_REQUIRED true)
