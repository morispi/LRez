cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(LRez CXX)
set(CMAKE_CXX_STANDARD 11)

find_package(Boost REQUIRED COMPONENTS iostreams)
find_package(ZLIB REQUIRED)
add_subdirectory(bamtools EXCLUDE_FROM_ALL)

file(GLOB SRC_FILES src/*.cpp src/subcommands/*.cpp)
add_executable(LRez ${SRC_FILES})
target_include_directories(LRez PRIVATE src/include ${BamTools_SOURCE_DIR}/src)
target_link_libraries(LRez BamTools Boost::iostreams ZLIB::ZLIB)
