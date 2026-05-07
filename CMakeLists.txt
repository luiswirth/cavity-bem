cmake_minimum_required(VERSION 3.25)
project(my_project LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# tell CMake where Homebrew's libomp lives
set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
set(OpenMP_CXX_LIB_NAMES "omp")
set(OpenMP_omp_LIBRARY /opt/homebrew/opt/libomp/lib/libomp.dylib)
include_directories(/opt/homebrew/opt/libomp/include)

find_package(OpenMP REQUIRED)

include(FetchContent)
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
)
FetchContent_MakeAvailable(eigen)
FetchContent_Declare(
  bembel
  GIT_REPOSITORY https://github.com/temf/bembel.git
  GIT_TAG v1.0
)
FetchContent_MakeAvailable(bembel)

add_executable(my_project src/main.cpp)
target_link_libraries(my_project PRIVATE Eigen3::Eigen OpenMP::OpenMP_CXX)
target_include_directories(my_project PRIVATE ${bembel_SOURCE_DIR})
