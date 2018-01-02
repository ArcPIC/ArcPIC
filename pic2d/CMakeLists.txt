# CMakeLists file for the core of ArcPIC
# Normally called by the CMakeLists file in the root of the repository

#Setup the target ArcPIC
file(GLOB ArcPIC_sources src/*.cpp)
add_executable(ArcPIC ${ArcPIC_sources})
target_include_directories(ArcPIC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/h)

#Use the SuperLU library
target_link_libraries(ArcPIC superlu)
message(${SuperLU_SOURCE_DIR})
target_include_directories(ArcPIC PUBLIC ${SuperLU_SOURCE_DIR}/SRC)

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
  target_include_directories(ArcPIC PUBLIC ${GSL_INCLUDE_DIRS})
  target_link_libraries(ArcPIC ${GSL_LIBRARIES})
else(OPENMP_FOUND)
  message("ERROR: OpenMP could not be found.")
endif(OPENMP_FOUND)
