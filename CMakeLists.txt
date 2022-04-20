# Fluid flow sim
#
# Copyright (c) 2022 Pascal Massimino (skal -> "pascal.massimino@gmail.com")
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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

# ##############################################################################
