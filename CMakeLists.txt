project(rsmd)
cmake_minimum_required(VERSION 3.6)


set( CMAKE_CXX_STANDARD 17 )

message( STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}" )
message( STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}" )
message( STATUS "BUILD TYPE: ${CMAKE_BUILD_TYPE}" )

if(APPLE)
    SET( CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "/usr/local/opt/boost@1.59/include" ) # required for finding boost on mac
    SET( CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "/usr/local/opt/boost@1.59/lib" )
endif()

find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})
  

# specific flags
set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g3 -ggdb -D_DEBUG -DDEBUG -frounding-math -fsignaling-nans")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g1 -ggdb -DNDEBUG -ffast-math")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -g0       -DNDEBUG -ffast-math")

# general flags
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -fmax-errors=3 -Weffc++")


# build executable from main
include_directories("./src/")   # in order to find includes
file( GLOB sources 
    src/*.cpp 
    src/*/*.cpp
)
message( STATUS "compiling: ${sources}")


add_executable( rsmd ${sources})


# link
target_link_libraries(rsmd ${STDCXX_LDFLAGS} "-lboost_program_options -lstdc++fs")

