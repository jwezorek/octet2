cmake_minimum_required(VERSION 3.24)

project(octet
    VERSION 0.1.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Boost 1.80 REQUIRED)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)

add_executable(octet
    src/main.cpp
    src/cell.cpp
    src/cell_set_polyhedron.cpp
    src/octet.cpp
    src/cell_set.cpp
    src/polyhedron.cpp
    src/color.cpp
    src/colored_cell_set_polyhedron.cpp
    src/colored_polyhedron.cpp
    src/colored_vector_dictionary.cpp
)

target_link_libraries(octet
    PRIVATE
        Boost::boost
        Eigen3::Eigen
)

target_compile_features(octet
    PRIVATE
        cxx_std_23
)

if(MSVC)
    target_compile_options(octet PRIVATE /W4 /permissive-)
else()
    target_compile_options(octet PRIVATE -Wall -Wextra -Wpedantic)
endif()