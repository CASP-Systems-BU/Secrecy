cmake_minimum_required(VERSION 3.19)

project(Secrecy)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-O3")

# Adding Dependencies
add_subdirectory(include/external-lib/sql-parser)

# Adding MPI Library
find_package(MPI REQUIRED)
if (MPI_FOUND)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
else (MPI_FOUND)
    message(SEND_ERROR "This application cannot compile without MPI")
endif (MPI_FOUND)

# Using libsodium
find_package(PkgConfig REQUIRED)
pkg_check_modules(_LIBSODIUM REQUIRED libsodium)

find_path(SODIUM_INCLUDE_DIR sodium.h HINTS ${LIBSODIUM_INCLUDE_DIRS} /usr/local/include /opt/local/include /opt/include)
find_library(SODIUM_LIBRARY NAMES sodium HINTS ${LIBSODIUM_LIBRARY_DIRS} /usr/local/lib /opt/local/lib /opt/lib)

# SecrecyLib
file(GLOB Secrecy_SOURCES src/core/*.cpp src/core/wrappers/*.cpp src/common/*.cpp src/benchmarking/*.cpp src/planner/*.cpp)
# add_library(SecrecyLib INTERFACE ${Secrecy_SOURCES})
# target_link_libraries(SecrecyLib INTERFACE ${MPI_LIBRARIES} ${SODIUM_LIBRARY})
# target_include_directories(SecrecyLib INTERFACE ${SODIUM_INCLUDE_DIR})

# TODO: factor out building all files with each executable.
# Adding testing executables
# src/tests/api/*.cpp 
file(GLOB TEST_EXECUTABLES src/tests/c_api/*.cpp src/tests/codegen/*.cpp)
foreach (_target ${TEST_EXECUTABLES})
    get_filename_component(fileName ${_target} NAME_WLE)
    add_executable(${fileName}.out ${_target} ${Secrecy_SOURCES})
    target_link_libraries(${fileName}.out PUBLIC ${MPI_LIBRARIES} ${SODIUM_LIBRARY} src)
    target_include_directories(${fileName}.out PUBLIC ${SODIUM_INCLUDE_DIR})
endforeach ()

# Adding experiments executables
file(GLOB TEST_EXECUTABLES src/experiments/*.cpp)
foreach (_target ${TEST_EXECUTABLES})
    get_filename_component(fileName ${_target} NAME_WLE)
    add_executable(${fileName}.out ${_target} ${Secrecy_SOURCES} include/planner/operators/table_function.h)
    target_link_libraries(${fileName}.out PUBLIC ${MPI_LIBRARIES} ${SODIUM_LIBRARY} src)
    target_include_directories(${fileName}.out PUBLIC ${SODIUM_INCLUDE_DIR})
endforeach ()

# Add Planner
#add_executable(planner.out src/tests/codegen/planner.cpp ${Secrecy_SOURCES} src/planner/dpTable.cpp src/planner/sql_expression.cpp src/planner/rule.cpp src/planner/codeGen.cpp)
#target_link_libraries(planner.out PUBLIC ${MPI_LIBRARIES} ${SODIUM_LIBRARY} src)
#target_include_directories(planner.out PUBLIC ${SODIUM_INCLUDE_DIR})