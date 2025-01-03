###########################################################################
## Authors: Sayak Kundu (sakundu@ucsd.edu), Zhiang Wang (zhw033@ucsd.edu)
##          Dooseok Yoon (d3yoon@ucsd.edu)
## Copyright (c) 2024, The Regents of the University of California
## All rights reserved.
##
## BSD 3-Clause License
##
## Copyright (c) 2019, James Cherry, Parallax Software, Inc.
## Copyright (c) 2018, SangGi Do and Mingyu Woo
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##
## * Redistributions of source code must retain the above copyright notice, this
##   list of conditions and the following disclaimer.
##
## * Redistributions in binary form must reproduce the above copyright notice,
##   this list of conditions and the following disclaimer in the documentation
##   and/or other materials provided with the distribution.
##
## * Neither the name of the copyright holder nor the names of its
##   contributors may be used to endorse or promote products derived from
##   this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
## ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
## LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
## CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
## SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
## INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
## CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
## ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
## POSSIBILITY OF SUCH DAMAGE.
###########################################################################

include("openroad")
find_package(TCL)

add_library(sa1d_lib
  src/Objects.cpp
  src/OptSA.cpp
  src/Worker.cpp
)

target_link_libraries(sa1d_lib
  PUBLIC
    odb
    OpenSTA
  PRIVATE
    utl_lib
)


target_include_directories(sa1d_lib
  PUBLIC
    include
)

# Add OpenMP support
find_package(OpenMP REQUIRED)
if (OpenMP_CXX_FOUND)
  target_link_libraries(sa1d_lib PRIVATE OpenMP::OpenMP_CXX)
endif()

swig_lib(NAME         sa1d
         NAMESPACE    sa1d
         I_FILE       src/OptSA.i
         SCRIPTS      src/OptSA.tcl
)

target_sources(sa1d
  PRIVATE
    src/MakeOptSA.cpp
)

target_include_directories(sa1d
  PUBLIC
    include
  PRIVATE
    src
)

target_link_libraries(sa1d
  PRIVATE
    ${TCL_LIBRARY}
    odb
    sa1d_lib
    OpenSTA
)

messages(
  TARGET sa1d
)

################################################################

if (Python3_FOUND AND BUILD_PYTHON)
  swig_lib(NAME          sa1d_py
           NAMESPACE     sa1d
           LANGUAGE      python
           I_FILE        src/OptSA-py.i
           SWIG_INCLUDES ${PROJECT_SOURCE_DIR}/include/sa1d
                         ${ODB_HOME}/src/swig/common
                         ${ODB_HOME}/src/swig/python
           SCRIPTS       ${CMAKE_CURRENT_BINARY_DIR}/sa1d_py.py
  )

  target_link_libraries(sa1d_py
    PUBLIC
      ${Boost_LIBRARIES}
      sa1d_lib
      odb
  )

endif()


