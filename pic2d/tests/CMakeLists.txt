# Do not omit !!!RELATIVE!!!
#file(GLOB_RECURSE pyfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
#create_symlinks(${pyfiles})

# COLLTEST is broken
# add_executable(colltest colltest.cpp)
# target_link_libraries(colltest ArcPICcore)
create_symlinks(colltestAna.py)

add_executable(rngtest rngtest.cpp)
target_link_libraries(rngtest ArcPICcore)
create_symlinks(rngtestAna.py)

add_executable(rngparatest rngparatest.cpp)
target_link_libraries(rngparatest ArcPICcore)

#add_executable(picFieldTest picFieldTest.cpp)
#target_link_librares(picFieldTest ArcPICcore)
create_symlinks(plotPicFieldTest.py)
