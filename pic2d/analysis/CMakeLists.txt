# Do not omit !!!RELATIVE!!!
file(GLOB_RECURSE pyfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
create_symlinks(${pyfiles})
