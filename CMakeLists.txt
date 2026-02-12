cmake_minimum_required(VERSION 3.16)
project(fast_patcher LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ── Compiler flags ─────────────────────────────────────────────────────────
# NOTE: Inside generator expressions, flags MUST be semicolon-separated
#       (spaces are treated as literal characters, not list separators).
set(cxx_flags
    $<$<CXX_COMPILER_ID:MSVC>:-D_SCL_SECURE_NO_WARNINGS;/openmp:experimental;/MP;/std:c++17;/bigobj>
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<CONFIG:Release>>:-Wall;-m64;-fopenmp;-O3;-Wno-unused-function>
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<CONFIG:Debug>>:-Wall;-m64;-fopenmp;-O0;-g3;-fno-inline;-Wno-unused-function>
    $<$<AND:$<CXX_COMPILER_ID:AppleClang>,$<CONFIG:Release>>:-Wall;-O3;-Wno-unused-function>
    $<$<AND:$<CXX_COMPILER_ID:AppleClang>,$<CONFIG:Debug>>:-Wall;-O0;-g3;-fno-inline;-fno-omit-frame-pointer;-Wno-unused-function>
    $<$<AND:$<CXX_COMPILER_ID:Clang>,$<NOT:$<CXX_COMPILER_ID:AppleClang>>,$<CONFIG:Release>>:-Wall;-m64;-fopenmp;-O3;-Wno-unused-function>
    $<$<AND:$<CXX_COMPILER_ID:Clang>,$<NOT:$<CXX_COMPILER_ID:AppleClang>>,$<CONFIG:Debug>>:-Wall;-m64;-fopenmp;-O0;-g3;-fno-inline;-Wno-unused-function>
)

add_library(fast_patcher_flags INTERFACE)
target_compile_options(fast_patcher_flags INTERFACE
    $<$<COMPILE_LANGUAGE:CXX>:${cxx_flags}>
)

# Ensure CMake's built-in Debug flags also use -O0 -g3
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3" CACHE STRING "C++ Debug flags" FORCE)
set(CMAKE_C_FLAGS_DEBUG "-O0 -g3" CACHE STRING "C Debug flags" FORCE)

# ── Third-party dependencies (fetched via cmake/recipes/) ──────────────────
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/recipes")

# Find OpenMP (required for the openmp example)
# macOS with AppleClang needs hints to find Homebrew's libomp
if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    execute_process(
        COMMAND brew --prefix libomp
        OUTPUT_VARIABLE LIBOMP_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(LIBOMP_PREFIX)
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY "${LIBOMP_PREFIX}/lib/libomp.dylib")
    endif()
endif()
find_package(OpenMP REQUIRED)

include(eigen)
include(libigl)
include(spdlog)
include(cli11)

# ── Library: fast_patcher_lib (Lloyd clustering + helpers) ─────────────────
file(GLOB LIB_SOURCES src/*.cpp)
add_library(fast_patcher_lib STATIC ${LIB_SOURCES})
target_include_directories(fast_patcher_lib PUBLIC src/)
target_link_libraries(fast_patcher_lib PUBLIC fast_patcher_flags spdlog::spdlog)

# Debug/Release definitions: DEBUG defined in Debug, NDEBUG in Release
target_compile_definitions(fast_patcher_lib PUBLIC
    $<$<CONFIG:Debug>:DEBUG>
    $<$<CONFIG:Release>:NDEBUG>
)

# ── Demo executable ────────────────────────────────────────────────────────
add_executable(test_cluster demo/test_cluster.cpp)
target_include_directories(test_cluster PRIVATE demo/utils/)
target_link_libraries(test_cluster PRIVATE
    fast_patcher_lib
    igl::core
    Eigen3::Eigen
    spdlog::spdlog
    CLI11::CLI11
)
