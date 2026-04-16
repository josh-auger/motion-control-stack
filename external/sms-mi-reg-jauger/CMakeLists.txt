PROJECT(SMSMIREG CXX)

cmake_minimum_required(VERSION 3.14)

# Environment requirements:
# OpenMP has to be installed
# ITK_DIR should be set to point to the Insight Toolkit e.g. version 5.2
# NLopt_DIR should be set to point to the NLopt library version 2.7.0
#    Use lower case opt in the name NLopt_DIR.
# LD_LIBRARY_PATH needs to point to ITK.

if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17) # Supported values are ``11``, ``14``, and ``17``.
endif()
if(NOT CMAKE_CXX_STANDARD_REQUIRED)
  set(CMAKE_CXX_STANDARD_REQUIRED YES)
endif()
if(NOT CMAKE_CXX_EXTENSIONS)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()
set(VALID_CXX_STANDARDS "11" "14" "17")
if(NOT CMAKE_CXX_STANDARD IN_LIST VALID_CXX_STANDARDS )
   MESSAGE(FATAL_ERROR "CMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD} not in known standards list\n ${VALID_CXX_STANDARDS} .")
endif()

set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_INSTALL_PREFIX "." )


FIND_PACKAGE(ITK REQUIRED )
INCLUDE(${ITK_USE_FILE})

LINK_LIBRARIES(${ITK_LIBRARIES})

# Add the quill package for logging.
# For system installation of quill:
find_package(quill REQUIRED)
# For local quill:
# add_subdirectory(quill)

ADD_EXECUTABLE(sms-mi-reg sms-mi-reg.cpp)

MESSAGE(STATUS "Make sure the compiler is appropriate for the platform. " )
# If you are running this on a Mac, set:
#set(CMAKE_C_COMPILER clang)
#set(CMAKE_CXX_COMPILER clang++)

# If you are running on linux, use devtoolset-8 and gcc 8 or newer
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

if(CMAKE_C_COMPILER_ID MATCHES "Clang\$")
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY omp)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang\$")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_CXX_INCLUDE_DIRS "/usr/local/include")
    set(OpenMP_omp_LIBRARY omp)
else()
MESSAGE("CMAKE_CXX_COMPILER_ID IS ${CMAKE_CXX_COMPILER_ID}")
endif()

find_package(OpenMP REQUIRED)

if (OpenMP_CXX_FOUND)
    link_directories("${OpenMP_CXX_LIBRARIES}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

# The search path for NLOPT can be modified by setting e.g. NLOPT_DIR
FIND_PACKAGE( NLopt  2.7 REQUIRED CONFIG
HINTS ${NLOPT_DIR}
)
MESSAGE( NOTICE "NLOPT_DIR is now set to [${NLOPT_DIR}]." )

MESSAGE( NOTICE "NLOPT_INCLUDE_DIRS is set to [${NLOPT_INCLUDE_DIRS}]." )
MESSAGE( NOTICE "NLOPT_LIBRARIES is set to [${NLOPT_LIBRARIES}]." )
MESSAGE( NOTICE "NLOPT_LIBRARY_DIRS is set to [${NLOPT_LIBRARY_DIRS}]." )


#
#  target_link_libraries(cache-line-sharing PUBLIC OpenMP::OpenMP_CXX)
#  target_link_directories(cache-line-sharing 
#        PUBLIC /usr/local/lib          )
#  target_include_directories(cache-line-sharing PUBLIC /usr/local/include
#    PUBLIC "${OpenMP_CXX_INCLUDE_DIRS}" )


  target_link_libraries(sms-mi-reg 
      PUBLIC OpenMP::OpenMP_CXX
      PUBLIC ${NLOPT_LIBRARIES}
      PRIVATE quill::quill
  )
  target_link_directories(sms-mi-reg
        PUBLIC /usr/local/lib 
        INTERFACE ${NLOPT_LIBRARY_DIRS} 
  )
  target_include_directories(sms-mi-reg PUBLIC /usr/local/include
    PUBLIC "${OpenMP_CXX_INCLUDE_DIRS}" 
    PUBLIC ${NLOPT_INCLUDE_DIRS} 
  )

#ADD_EXECUTABLE(volume-ply-reg volume-ply-reg.cpp)
#target_link_libraries(volume-ply-reg 
#      PUBLIC OpenMP::OpenMP_CXX
#      PUBLIC ${NLOPT_LIBRARIES}
#)
#target_link_directories(volume-ply-reg
#        PUBLIC /usr/local/lib 
#        INTERFACE ${NLOPT_LIBRARY_DIRS} 
#)
#target_include_directories(volume-ply-reg PUBLIC /usr/local/include
#    PUBLIC "${OpenMP_CXX_INCLUDE_DIRS}" 
#    PUBLIC ${NLOPT_INCLUDE_DIRS} 
#)

ADD_EXECUTABLE(image-to-ply-converter image-to-ply-converter.cpp)
target_link_libraries(image-to-ply-converter ${ITK_LIBRARIES})


endif(OpenMP_CXX_FOUND)

ADD_EXECUTABLE(polar-decomposition polar-decomposition.cpp)
ADD_EXECUTABLE(writeAffineTransform writeAffineTransform.cpp)

target_link_libraries(writeAffineTransform
      PUBLIC ${ITK_LIBRARIES}
  )

target_include_directories(writeAffineTransform 
    PUBLIC /usr/local/include
    PUBLIC ${ITK_INCLUDE_DIRS} 
)


