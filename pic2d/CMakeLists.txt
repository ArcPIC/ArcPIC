# CMakeLists file for the core of ArcPIC
# Normally called by the CMakeLists file in the root of the repository

include(${CMAKE_CURRENT_SOURCE_DIR}/create_symlinks.cmake)

#Setup the target ArcPIC
#file(GLOB ArcPIC_sources src/*.cpp)
#add_executable(ArcPIC ${ArcPIC_sources})
#target_include_directories(ArcPIC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/h)

#Setup a ArcPIC core library
file(GLOB ArcPIC_sources_nomain src/*.cpp)
list(REMOVE_ITEM ArcPIC_sources_nomain ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
add_library(ArcPICcore STATIC ${ArcPIC_sources_nomain})
target_include_directories(ArcPICcore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/h)

## Link external libraries ##

#Use the SuperLU library
target_link_libraries(ArcPICcore superlu)
target_include_directories(ArcPICcore PUBLIC ${SuperLU_SOURCE_DIR}/SRC)

#Use OpenMP for paralellization
include(FindOpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else(OPENMP_FOUND)
  message("ERROR: OpenMP could not be found.")
endif(OPENMP_FOUND)

#Use GSL for RNGs
include(FindGSL)
if(GSL_FOUND)
  target_include_directories(ArcPICcore PUBLIC ${GSL_INCLUDE_DIRS})
  target_link_libraries(ArcPICcore ${GSL_LIBRARIES})
else(GSL_FOUND)
  message("ERROR: LibGSL could not be found.")
endif(GSL_FOUND)

#Use HDF5 for the heavy IO
find_package (HDF5 COMPONENTS CXX) # Find non-cmake built HDF5
target_include_directories (ArcPICcore PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(ArcPICcore ${HDF5_CXX_LIBRARIES})

## Executables

#Main executable
add_executable(ArcPIC src/main.cpp)
target_link_libraries(ArcPIC ArcPICcore)

## Other stuff ##

#Copy the input file
configure_file(input.txt input.txt COPYONLY)

#Get the analyses
add_subdirectory(analysis)
#Get the tests
add_subdirectory(tests)
