cmake_minimum_required(VERSION 3.10.0)
project (CATViz)

#find_package(Boost REQUIRED)

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/deps/glm
    )

link_directories(
    )

add_executable(cviz src/main.cpp)
if (MSVC)
    set_directory_properties(PROPERTIES 
        COMPILE_DEFINITIONS _VARIADIC_MAX=10
        ) 
endif()
