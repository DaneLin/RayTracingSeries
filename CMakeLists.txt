cmake_minimum_required(VERSION 3.0.0)
project(PBR VERSION 0.1.0 LANGUAGES C CXX)

# FIND OPENMP
FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)
    message(STATUS "found openmp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    #set(CMAKE_CXX_FLAGS"${CMAKE_CXX_FLAGS} -fopenmp")
    message(STATUS "optional:-fopenmp")
endif()

set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}")

file(COPY ${SOURCE_DIR}/images DESTINATION ${BINARY_DIR})

include_directories(${PROJECT_PATH}includes)
include_directories(${PROJECT_PATH}external)
include_directories(${PROJECT_PATH}images)


add_executable(PBR main.cpp)

target_link_libraries(PBR ${OpenMP_LIBS})

