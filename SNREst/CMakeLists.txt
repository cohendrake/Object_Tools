cmake_minimum_required(VERSION 2.8)

project(SNREst)

#project path

set(root_dir ".")

#set debug & release build option

set(CMAKE_C_FLAGS_DEBUG "$ENV{CFLAGS} -DINFO_OUT -O0 -Wall -g -ggdb")

set(CMAKE_C_FLAGS_RELEASE "$ENV{CFLAGS} -O3 -Wall")

#build executable

set(target_out "SNREst_${CMAKE_BUILD_TYPE}")

include_directories(${root_dir}/Include)

file(GLOB src_c ${root_dir}/Src/*.c)

add_executable(${root_dir}/Bin/${target_out} ${src_c})

TARGET_LINK_LIBRARIES(${root_dir}/Bin/${target_out} m) #link to math libarary

