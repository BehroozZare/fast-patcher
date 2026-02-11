cmake_minimum_required(VERSION 3.16)
project(fast_patcher LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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
