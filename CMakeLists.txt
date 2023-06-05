cmake_minimum_required(VERSION 3.1)

#set(CMAKE_C_COMPILER "/unicluster/opt/gcc/8/bin/gcc")
#set(CMAKE_CXX_COMPILER "/unicluster/opt/gcc/8/bin/g++")

project(compressibleEulerDGMFEM)

### include sources
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

### CXX flags
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -g")

### include hypre & metis libs for parallel computing
set(HYPRE_DIR "/home/vkorchagova/software/hypre/src/hypre")
set(METIS_DIR "/usr/include")

message(STATUS "HYPRE dir: ${HYPRE_DIR}")
message(STATUS "METIS dir: ${METIS_DIR}")
find_library(LIB_METIS REQUIRED NAMES METIS metis metis.a HINTS "${METIS_DIR}/../lib")
message(STATUS "METIS lib file: ${LIB_METIS}")


### include MFEM library
message(STATUS "Looking for mfem ...")
set(MFEM_SOURCE_DIR $ENV{HOME}/software/mfem/build/usr/local/)
set(mfem_DIR ${MFEM_SOURCE_DIR}/lib/cmake/mfem)
set (LIB_MFEM_DIR ${mfem_DIR})
find_library(LIB_MFEM REQUIRED NAMES MFEM mfem mfem.a HINTS "${MFEM_SOURCE_DIR}/lib")
message(STATUS "Found mfem config in: ${mfem_DIR} (version ${MFEM_VERSION})")
message(STATUS "Lib file: ${LIB_MFEM}")
include_directories(${MFEM_SOURCE_DIR}/include ${HYPRE_DIR}/include ${dirlist})

### include libprofiler to profile the application
#find_library(LIB_PROFILER NAMES profiler libprofiler.so libprofiler.so.0.4.9 HINTS "/unicluster/home/viktoria.korchagova/software/gperftools/lib")
#message(STATUS "gperftools lib file: ${LIB_PROFILER}")
#find_library(LIB_TCMALLOC NAMES tcmalloc libtcmalloc.so HINTS "/unicluster/home/viktoria.korchagova/software/gperftools/lib")
#message(STATUS "TCMALLOC lib file: ${LIB_TCMALLOC}")


### find MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
message(STATUS "MPI libs: ${MPI_CXX_LIBRARIES} ")
message(STATUS "MPI incs: ${MPI_INCLUDE_PATH} ")
### find yaml lib for case interface
include_directories(/usr/include)
include_directories(/home/vkorchagova/software/rapidyaml/build/usr/local/include)
include_directories(/usr/lib/cmake)
#find_package(LIBRYML REQUIRED)
#message(STATUS "ryml libs: ${LIB_RYML} ")


### complete linking
add_executable(${PROJECT_NAME} ${ALLFILES})
#target_link_libraries(${PROJECT_NAME} /usr/lib/libryml.so ${LIB_MFEM} -L${HYPRE_DIR}/../lib -lHYPRE -lrt ${MPI_LIBRARIES} ${LIB_METIS} -lstdc++fs)
target_link_libraries(${PROJECT_NAME} /home/vkorchagova/software/rapidyaml/build/usr/local/lib/libryml.a ${LIB_MFEM} -L${HYPRE_DIR}/../lib -lHYPRE -lrt ${MPI_LIBRARIES} ${LIB_METIS} -lstdc++fs)
#target_link_libraries(${PROJECT_NAME} $ENV{HOME}/software/rapidyaml-0.1.0/build/lib/libryml.a ${LIB_MFEM} -L${HYPRE_DIR}/../lib -lHYPRE   -lrt ${MPI_LIBRARIES} ${LIB_METIS} -lstdc++fs)


