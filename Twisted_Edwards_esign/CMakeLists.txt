cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

set(NAMESPACE Twisted_Edwards)
set(COMPONENT esign)
project(${NAMESPACE}_${COMPONENT}
    VERSION 0.0.0
    LANGUAGES C
)

add_library(${PROJECT_NAME}
    esign.h
    esign.c
)
find_package(Twisted_Edwards_ecurve 0.0.0 REQUIRED)
target_link_libraries(${PROJECT_NAME} Twisted_Edwards::ecurve)
add_library(${NAMESPACE}::${COMPONENT} ALIAS ${PROJECT_NAME})
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")

