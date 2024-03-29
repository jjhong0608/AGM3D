cmake_minimum_required(VERSION 3.18)
project(AGM3D)

set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

set(MKL_DIR "/opt/intel/compilers_and_libraries_2020.3.275/linux/mkl/include")
include_directories(${MKL_DIR} Header Header/Eigen)

add_executable(AGM3D main.cpp Header/util.cpp Header/util.h Header/Greenfunction.cpp Header/Greenfunction.h Header/GreenfunctionReactionDiffusion.cpp Header/GreenfunctionReactionDiffusion.h Header/coordinate.cpp Header/coordinate.h Header/value.cpp Header/value.h Header/matrixRow.cpp Header/matrixRow.h Header/axialLine.cpp Header/axialLine.h Header/plane.cpp Header/plane.h Header/readFile.cpp Header/readFile.h Header/point.cpp Header/point.h Header/solver.cpp Header/solver.h Header/matrix.cpp Header/matrix.h Header/writeFile.cpp Header/writeFile.h Header/pointHeat.cpp Header/pointHeat.h Header/function.cpp Header/function.h Header/matrixMulti.cpp Header/matrixMulti.h Header/matrixNormal.cpp Header/matrixNormal.h Header/writeFileMultiple.cpp Header/writeFileMultiple.h)

find_library(
        ILP64_MODULE
        NAMES libmkl_intel_lp64.so
        HINTS /opt/intel/compilers_and_libraries_2020.3.275/linux/mkl/lib/intel64
        REQUIRED
)

find_library(
        THREAD_MODULE
        NAMES libmkl_intel_thread.so
        HINTS /opt/intel/compilers_and_libraries_2020.3.275/linux/mkl/lib/intel64
        REQUIRED
)

find_library(
        CORE_MODULE
        NAMES libmkl_core.so
        HINTS /opt/intel/compilers_and_libraries_2020.3.275/linux/mkl/lib/intel64
        REQUIRED
)

find_library(
        SEQUENTIAL_MODULE
        NAMES libmkl_sequential.so
        HINTS /opt/intel/compilers_and_libraries_2020.3.275/linux/mkl/lib/intel64
        REQUIRED
)

find_library(
        IOMP_MODULE
        NAMES libiomp5.a
        HINTS /opt/intel/compilers_and_libraries_2020.3.275/linux/compiler/lib/intel64
        REQUIRED
)

find_library(
        LPTHREAD
        NAMES pthread
)

find_library(
        LM
        NAMES m
)

find_library(
        LDL
        NAMES dl
)

target_link_libraries(AGM3D ${ILP64_MODULE} ${THREAD_MODULE} ${CORE_MODULE} ${SEQUENTIAL_MODULE} ${IOMP_MODULE} ${LPTHREAD} ${LM} ${LDL})

target_compile_options(AGM3D PUBLIC -xCOMMON-AVX512)