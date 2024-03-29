cmake_minimum_required(VERSION 3.22.1)

project(boyle_pnc
  LANGUAGES
    C
    CXX
)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_C_FLAGS_RELEASE "-O2")
set(CMAKE_C_FLAGS_DEBUG "-Og")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS_DEBUG "-Og")

set(CMAKE_COMPILE_WARNING_AS_ERROR True)

add_compile_options(
    -Wall
    -Wextra
    -Wpedantic
    -Wno-uninitialized
    -Wno-unused-parameter
)

add_link_options(
    -fno-ipa-icf
)

set(CPM_SOURCE_CACHE third_party)
set(CPM_USE_LOCAL_PACKAGES True)

add_compile_definitions(
    SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE
    DOCTEST_CONFIG_SUPER_FAST_ASSERTS
)

include(CTest)
enable_testing()

include(cmake/CPM.cmake)
include(cmake/third_party.cmake)
include(cmake/utils.cmake)

include_directories(${PROJECT_SOURCE_DIR})

add_subdirectory(common)
add_subdirectory(kinetics)
add_subdirectory(math)
add_subdirectory(tasks)
add_subdirectory(workspaces)
