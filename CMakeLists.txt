cmake_minimum_required(VERSION 3.21)

project(main)

set(CMAKE_CXX_STANDARD 23)

add_executable(main main.cpp)

add_subdirectory(src)

target_link_libraries(main lp)

