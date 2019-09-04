cmake_minimum_required(VERSION 3.10.0)
project (CATViz)

find_package(Boost REQUIRED)

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories(
    ${Boost_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/deps/svgpp/include
    )

link_directories(
    ${Boost_LIBRARY_DIRS}
    )

add_executable(cviz src/main.cpp)
if (MSVC)
    set_directory_properties(PROPERTIES 
        COMPILE_DEFINITIONS _VARIADIC_MAX=10
        ) 
endif()
