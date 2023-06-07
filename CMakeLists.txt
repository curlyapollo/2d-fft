cmake_minimum_required(VERSION 3.14)

project(2d-fft)

set(CMAKE_CXX_STANDARD 20)
add_compile_options(-Wall -Werror -O2 -ffast-math)


Include(FetchContent)
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.3.2
)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.13.0
)
FetchContent_Declare(
        googlebenchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG v1.8.0
)
FetchContent_MakeAvailable(
        Catch2
        googletest
        googlebenchmark)

add_executable(tests test.cpp src/BitReverse.cpp)
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)
target_include_directories(tests PRIVATE src/include)

add_executable(bench bench.cpp src/BitReverse.cpp)
target_link_libraries(bench PRIVATE benchmark::benchmark)
target_include_directories(bench PRIVATE src/include)
