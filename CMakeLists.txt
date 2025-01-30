cmake_minimum_required(VERSION 3.17)
project(09-mc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

file(GLOB_RECURSE SRC_FILES
    "${CMAKE_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_SOURCE_DIR}/src/*/.h"
)
add_executable(09-mc ${SRC_FILES})

include(FetchContent)

# Download Eigen
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
)
FetchContent_MakeAvailable(eigen)

# Download Qhull
FetchContent_Declare(
    qhull
    GIT_REPOSITORY https://github.com/qhull/qhull.git
    GIT_TAG master
)
FetchContent_MakeAvailable(qhull)

# Find Threads
find_package(Threads REQUIRED)

# Link libraries
target_link_libraries(09-mc PRIVATE
    Threads::Threads
    Eigen3::Eigen
    qhullstatic_r
    qhullcpp
)

# Include directories
target_include_directories(09-mc PRIVATE
    "${CMAKE_SOURCE_DIR}/src"
    ${qhull_SOURCE_DIR}/src
)

# Optimized compiler options
target_compile_options(09-mc PRIVATE -Ofast -march=native)
