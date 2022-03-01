cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(lima LANGUAGES CXX CUDA)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 20)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++20")
set(source_dir "${PROJECT_SOURCE_DIR}/src/")



find_package(CUDA REQUIRED)
find_package(Chrono)


add_subdirectory(${PROJECT_SOURCE_DIR}/Dependencies/GLFW/)
find_package(OpenGL REQUIRED)

file(GLOB source_files "${source_dir}/*.c*")

add_executable(lima ${source_files})

target_link_libraries(lima glfw)
target_link_libraries(lima OpenGL::GL)

target_compile_features(lima PUBLIC cxx_std_20)
#set_target_properties(lima PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(lima PROPERTIES CUDA_ARCHITECTURES "61;75")


target_include_directories(lima PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(lima PUBLIC "${CUDA_INCLUDE_DIRS}")
