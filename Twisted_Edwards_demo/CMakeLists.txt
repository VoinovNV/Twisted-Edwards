cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

set(NAMESPACE Twisted_Edwards)
set(COMPONENT demo)
project(${NAMESPACE}_${COMPONENT}
    VERSION 0.0.0
    LANGUAGES C
)

add_executable(${PROJECT_NAME}
    Source.c
)

add_executable(${NAMESPACE}::${COMPONENT} ALIAS ${PROJECT_NAME})

find_package(Twisted_Edwards_ecurve 0.0.0 REQUIRED)
target_link_libraries(${PROJECT_NAME} Twisted_Edwards::ecurve)

find_package(Twisted_Edwards_esign 0.0.0 REQUIRED)
target_link_libraries(${PROJECT_NAME} Twisted_Edwards::esign)

find_package(Twisted_Edwards_tests 0.0.0 REQUIRED)
target_link_libraries(${PROJECT_NAME} Twisted_Edwards::tests)

find_library(LIBAKR libakrypt-static.a
    PATHS /home/user/Projects/Twisted_Edwards/libakrypt-0.x/
)
target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../libakrypt-0.x/source/>)
target_link_libraries(${PROJECT_NAME} ${LIBAKR})
