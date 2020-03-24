# Install script for directory: /home/user/Projects/Twisted_Edwards/libakrypt-0.x

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/user/Projects/Twisted_Edwards/build/libakrypt-static.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/libakrypt.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/ak_tools.h"
    "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/ak_buffer.h"
    "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/ak_oid.h"
    "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/ak_random.h"
    "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/ak_mpzn.h"
    "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/ak_curves.h"
    "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/ak_parameters.h"
    "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/ak_gf2n.h"
    "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/ak_network.h"
    "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/ak_hash.h"
    "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/ak_skey.h"
    "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/ak_hmac.h"
    "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/ak_mac.h"
    "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/ak_context_manager.h"
    "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/ak_sign.h"
    "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/ak_bckey.h"
    "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/ak_omac.h"
    "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/ak_mgm.h"
    "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/ak_key_manager.h"
    "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/asn_processor/ak_asn_codec_new.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/etc/libakrypt.conf")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/etc" TYPE FILE FILES "/home/user/Projects/Twisted_Edwards/libakrypt-0.x/libakrypt.conf")
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/user/Projects/Twisted_Edwards/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
