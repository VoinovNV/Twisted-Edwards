cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(Twisted_Edwards
    VERSION 0.0.0
    LANGUAGES CXX
)

set(SUBPROJECTS
    Twisted_Edwards_ecurve
    Twisted_Edwards_esign
    Twisted_Edwards_tests
    Twisted_Edwards_demo
)

macro(find_package)
    if(NOT "${ARGV0}" IN_LIST SUBPROJECTS)
        _find_package(${ARGV})
    endif()
endmacro()

foreach(subdir IN LISTS SUBPROJECTS)
    add_subdirectory(${subdir})
endforeach()
