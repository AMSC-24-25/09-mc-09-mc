cmake_minimum_required(VERSION 3.17)
project(09-mc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

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
    GIT_TAG master  # Usa 'master' o un tag di versione valido
)

set(QHULL_BUILD_SHARED OFF CACHE BOOL "" FORCE)
set(QHULL_BUILD_STATIC ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(qhull)

file(GLOB_RECURSE SRC_FILES
    "${CMAKE_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_SOURCE_DIR}/src/*/.cpp"
    "${CMAKE_SOURCE_DIR}/src/*/.h"
)

add_executable(09-mc ${SRC_FILES})
add_dependencies(09-mc qhull)

# Link package pthread
find_package(Threads REQUIRED)
target_link_libraries(09-mc PRIVATE Threads::Threads)

# Include directories
target_include_directories(09-mc PRIVATE
    "${CMAKE_SOURCE_DIR}/src"  # Include le directory del tuo progetto
    ${Eigen3_INCLUDE_DIRS}
    ${qhull_SOURCE_DIR}/src/libqhull_r
    ${qhull_SOURCE_DIR}/src/libqhullcpp
    ${qhull_SOURCE_DIR}/src
)

# Link libs
target_link_libraries(09-mc PRIVATE Eigen3::Eigen
   	qhullstatic_r  # libreria C “robusta” di Qhull
    	qhullcpp       # libreria C++ di Qhull
	
)

# Optimized compiler options
target_compile_options(09-mc PRIVATE -Ofast -march=native)