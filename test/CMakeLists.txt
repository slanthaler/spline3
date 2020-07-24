cmake_minimum_required(VERSION 3.10)

file(GLOB_RECURSE sources "src/*.cpp")
file(GLOB_RECURSE headers "src/*.hpp" "include/*.hpp")

set(link_dir "/home/sam/Progs/play/c++/spline3/include")

add_executable(TestAll main_TestAll.cpp ${sources} ${headers})
include_directories(TestAll
  PRIVATE
  link_dir
  ${CMAKE_CURRENT_SOURCE_DIR}/include/)
target_link_libraries(TestAll spline3 GTest::GTest GTest::Main)
