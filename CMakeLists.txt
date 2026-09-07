cmake_minimum_required(VERSION 3.20)
project(numsolve CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(nlohmann_json CONFIG REQUIRED)

add_library(numsolve_core
    src/matrix.cpp
    src/elimination.cpp
    src/linalg.cpp
    src/condition.cpp
    src/diagnostics.cpp
    src/parser.cpp
    src/report.cpp
)
target_include_directories(numsolve_core PUBLIC include)
target_link_libraries(numsolve_core PUBLIC nlohmann_json::nlohmann_json)
if(MINGW)
    target_link_options(numsolve_core PUBLIC -static)
endif()

add_executable(numsolve src/main.cpp)
target_link_libraries(numsolve PRIVATE numsolve_core)

option(NUMSOLVE_BUILD_TESTS "Build the numsolve unit test suite" ON)
if(NUMSOLVE_BUILD_TESTS)
    find_package(Catch2 3 CONFIG REQUIRED)
    add_executable(numsolve_tests
        tests/test_matrix.cpp
        tests/test_elimination.cpp
        tests/test_condition.cpp
        tests/test_parser.cpp
    )
    target_link_libraries(numsolve_tests PRIVATE numsolve_core Catch2::Catch2WithMain)

    enable_testing()
    include(CTest)
    add_test(NAME numsolve_tests COMMAND numsolve_tests)
endif()
