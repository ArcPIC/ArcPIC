# support/CMakeLists.txt
# Collects the various libraries that we don't expect to find installed on the system

#Make sure the submodules are all there
include(submodules.cmake)

#SuperLU library
set(enable_tests CACHE BOOL OFF)
add_subdirectory(superlu)
