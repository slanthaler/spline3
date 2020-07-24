###
cmake_minimum_required(VERSION 3.10)

### Using Hunter package manager ###
include("cmake/HunterGate.cmake")
HunterGate(
  URL "https://github.com/ruslo/hunter/archive/v0.23.111.tar.gz"
  SHA1 "4189bbd869abe4eb69b5303d3f916c03d34dc154"
)

### define project ###
project(CmakeUnitTests
  VERSION 0.1
  DESCRIPTION "This is a example project for using CMake and the GoogleTest framework"
  LANGUAGES CXX)

### Download dependencies ###
hunter_add_package(GTest)

### Find dependencies ###
find_package(GTest REQUIRED)


### project details ###
# if this is the main project
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  # ensure -std=c++xx instead of -std=g++xx
  set(CMAKE_CXX_EXTENSIONS OFF) 

  # Docs, if this is the main app
  find_package(Doxygen)
  if(Doxygen_FOUND)
    add_subdirectory(docs)
  else()
    message(STATUS "Doxygen not found, not building docs")
  endif()
  
endif()

### Require out-of-source builds ###
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

### create a library with all routines ###
file(GLOB_RECURSE sources "${CMAKE_SOURCE_DIR}/src/*.cpp")
file(GLOB_RECURSE headers "${CMAKE_SOURCE_DIR}/include/*.hpp")

## following for library which needs to be compiled!
#add_library(spline3 STATIC ${sources} ${headers})
#target_include_directories(spline3 PUBLIC ${CMAKE_SOURCE_DIR}/include)
#set_target_properties(spline3 PROPERTIES ARCHIVE_OUTPUT_DIRECTORY lib)

## library header-only
add_library(spline3 INTERFACE)
target_include_directories(spline3 INTERFACE ${CMAKE_SOURCE_DIR}/include)

##
#add_executable(quick_test test.cpp ${sources} ${headers})
#target_include_directories(quick_test PRIVATE ${CMAKE_SOURCE_DIR}/include)


### testing (googletest) ###
enable_testing()
add_subdirectory(test)



  
