cmake_minimum_required(VERSION 3.1)

#set(CMAKE_C_COMPILER "/unicluster/opt/gcc/8/bin/gcc")
#set(CMAKE_CXX_COMPILER "/unicluster/opt/gcc/8/bin/g++")

project(dgns)


FILE(GLOB ALLFILES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.hpp
  )

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${curdir}/${child})
    ENDIF()
  ENDFOREACH()
  LIST(APPEND dirlist ${curdir}/${child})
  SET(${result} ${dirlist})
ENDMACRO()


SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR}/src)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2")

# METIS
set(METIS_DIR "$ENV{HOME}/software/metis-5.1.0")
set(METIS_LIB_DIR "${METIS_DIR}/lib")
find_library(LIB_METIS REQUIRED NAMES METIS metis metis.a libmetis libmetis.a HINTS "${METIS_LIB_DIR}")
message(STATUS "Found METIS config in: ${METIS_LIB_DIR}")
message(STATUS "Lib file: ${LIB_METIS}")

# HYPRE
set(HYPRE_DIR "$ENV{HOME}/software/hypre/src/hypre/include")
set(HYPRE_LIB_DIR "$ENV{HOME}/software/hypre/src/hypre/lib")
find_library(LIB_HYPRE REQUIRED NAMES HYPRE HYPRE.a libHYPRE.a HINTS "${HYPRE_LIB_DIR}")
message(STATUS "Found HYPRE config in: ${HYPRE_LIB_DIR}")
message(STATUS "Lib file: ${LIB_HYPRE}")

# MFEM
message(STATUS "Looking for mfem ...")
set(MFEM_SOURCE_DIR $ENV{MFEM_ROOT})
set(mfem_DIR $ENV{MFEM_ROOT}/lib/cmake/mfem)
set (LIB_MFEM_DIR ${mfem_DIR})
find_library(LIB_MFEM REQUIRED NAMES MFEM mfem mfem.a HINTS "~/software/mfem/build/usr/local/lib")
message(STATUS "Found mfem config in: ${mfem_DIR} (version ${MFEM_VERSION})")
message(STATUS "Lib file: ${LIB_MFEM}")


#find_library(LIB_PROFILER NAMES profiler libprofiler.so libprofiler.so.0.4.9 HINTS "/unicluster/home/viktoria.korchagova/extendedSrc/gperftools/lib")

#message(STATUS "gperftools lib file: ${LIB_PROFILER}")



include_directories(${MFEM_SOURCE_DIR}/include ${HYPRE_DIR} ${dirlist})

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

message(STATUS "MPI libs: ${MPI_CXX_LIBRARIES} ")
message(STATUS "MPI incs: ${MPI_INCLUDE_PATH} ")

add_executable(${PROJECT_NAME} ${ALLFILES})

target_link_libraries(${PROJECT_NAME} ${LIB_MFEM} ${LIB_HYPRE} ${LIB_METIS} -lrt ${MPI_LIBRARIES})