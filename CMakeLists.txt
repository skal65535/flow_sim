# Fluid flow sim

cmake_minimum_required(VERSION 3.7)
project(
  flow_sim
  LANGUAGES C CXX
  VERSION 0.1.0)
set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Build type: Release, Debug or RelWithDebInfo" FORCE)
endif()

# ##############################################################################
# target

find_package(SDL)

add_executable(sim ${CMAKE_CURRENT_SOURCE_DIR}/src/sim.cc)
if (SDL_FOUND)
  target_link_libraries(sim ${SDL_LIBRARY})
  target_include_directories(sim PRIVATE ${SDL_INCLUDE_DIR})
  target_compile_definitions(sim PUBLIC SIM_HAVE_SDL)
  set(SIM_HAVE_SDL 1)
endif()


# small tool to compute filters
add_executable(gauss ${CMAKE_CURRENT_SOURCE_DIR}/src/gauss.cc)
