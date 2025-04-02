cmake_minimum_required(VERSION 3.25)

project(
  Knots
  VERSION 1.0
  LANGUAGES CXX)

set(CMAKE_INSTALL_PREFIX
    "${CMAKE_BINARY_DIR}/install"
    CACHE PATH "CMAKE_INSTALL_PREFIX overwrite as requested by Abseil" FORCE)

# A workaround since MacOSX15 doesn't seem to be compatible with GCC 14.2:
set(CMAKE_OSX_SYSROOT /Library/Developer/CommandLineTools/SDKs/MacOSX14.sdk)

# ##############################################################################
# Compile options                                                              #
# ##############################################################################
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")
if(NOT gcc_like_cxx)
  message(
    FATAL_ERROR
      "Unsupported complier; please specify its flags in the CMakeLists.txt file."
  )
endif()

add_library(compiler_flags INTERFACE)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  message("Detected debug mode.")
  target_compile_options(
    compiler_flags
    INTERFACE
      "-Wall;-Wextra;-Wshadow;-Wformat=2;-Wunused;-Wno-sign-compare")

elseif(CMAKE_BUILD_TYPE STREQUAL Release)
  message("Detected release mode.")
elseif(CMAKE_BUILD_TYPE STREQUAL RelWithDebInfo)
  message("Detected release mode with debug info.")
  # You may want to include -fno-inline to help pinpointing the offending
  # function:
  #
  # target_compile_options(compiler_flags INTERFACE "-fno-inline")
else()
  message(
    FATAL_ERROR
      "Unknown build type (${CMAKE_BUILD_TYPE}); please add its configuration to the CMakeLists.txt file."
  )
endif()

# ##############################################################################
# Abseil                                                                       #
# ##############################################################################
include(FetchContent)
set(ABSL_ENABLE_INSTALL ON)
set(ABSL_PROPAGATE_CXX_STD ON)
FetchContent_Declare(
  absl
  GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git
  GIT_TAG origin/master
  OVERRIDE_FIND_PACKAGE)
FetchContent_MakeAvailable(absl)

# ##############################################################################
# Libraries                                                                    #
# ##############################################################################
add_library(sparse_bool_matrix_utils STATIC sparse_bool_matrix_utils.cpp)
target_link_libraries(
  sparse_bool_matrix_utils
  PUBLIC compiler_flags
  PRIVATE sparse_bool_matrix absl::check)

add_library(sparse_bool_matrix STATIC sparse_bool_matrix.cpp)
target_link_libraries(
  sparse_bool_matrix
  PUBLIC compiler_flags absl::flat_hash_set
  PRIVATE sparse_bool_matrix_utils sparse_low_memory_bool_matrix)

add_library(sparse_low_memory_bool_matrix STATIC
            sparse_low_memory_bool_matrix.cpp)
target_link_libraries(
  sparse_low_memory_bool_matrix
  PUBLIC compiler_flags sparse_bool_matrix
  PRIVATE sparse_bool_matrix_utils absl::check)

add_library(sparse_bool_matrix_fast_column_operations STATIC
            sparse_bool_matrix_fast_column_operations.cpp)
target_link_libraries(
  sparse_bool_matrix_fast_column_operations
  PUBLIC compiler_flags sparse_bool_matrix
  PRIVATE sparse_bool_matrix_utils)

add_library(algebra STATIC algebra.cpp)
target_link_libraries(
  algebra
  PUBLIC compiler_flags sparse_bool_matrix
  PRIVATE sparse_bool_matrix_fast_column_operations
          sparse_low_memory_bool_matrix)

add_library(pdcode STATIC pdcode.cpp)
target_link_libraries(pdcode PUBLIC compiler_flags)

add_library(symmetry STATIC symmetry.cpp)
target_link_libraries(symmetry PUBLIC compiler_flags algebra chaincomplex
                                      pdcode)

add_library(chaincomplex STATIC chaincomplex.cpp)
target_link_libraries(chaincomplex PUBLIC compiler_flags algebra pdcode)

add_library(khovanov STATIC khovanov.cpp)
target_link_libraries(khovanov PUBLIC compiler_flags algebra chaincomplex
                                      pdcode symmetry)

# ##############################################################################
# Binaries                                                                     #
# ##############################################################################
add_executable(full_khovanov full_khovanov.cpp)
target_link_libraries(full_khovanov PRIVATE compiler_flags khovanov pdcode)

add_executable(action action.cpp)
target_link_libraries(action PRIVATE compiler_flags khovanov pdcode)

# ##############################################################################
# Tests                                                                        #
# ##############################################################################
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG origin/main)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt
    ON
    CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

include(CTest)
enable_testing()

add_library(test_util INTERFACE)
target_link_libraries(test_util INTERFACE GTest::gtest GTest::gmock)

add_executable(sparse_bool_matrix_test sparse_bool_matrix_test.cpp)
target_link_libraries(
  sparse_bool_matrix_test PRIVATE GTest::gtest_main GTest::gmock test_util
                                  sparse_bool_matrix)

add_executable(algebra_test algebra_test.cpp)
target_link_libraries(algebra_test PRIVATE GTest::gtest_main GTest::gmock
                                           algebra test_util)

include(GoogleTest)
gtest_discover_tests(sparse_bool_matrix_test)
gtest_discover_tests(algebra_test)
