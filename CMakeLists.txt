cmake_minimum_required(VERSION 3.20)
project(TB_2025_Oct VERSION 0.3.0 LANGUAGES CXX)

include(FetchContent)

# ---- C++ standard (ROOT on most distros is built with C++17) ----
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---- Build type ----
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ---- macOS-only settings (guarded) ----
if(APPLE)
  # Set to your macOS version; harmless on Arch since guarded
  set(CMAKE_OSX_DEPLOYMENT_TARGET "14.0" CACHE STRING "" FORCE)

  # AppleClang already uses libc++; don't force flags here.
  # Optional: ensure AppleClang when on macOS
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    message(WARNING "On macOS, prefer /usr/bin/clang++ (AppleClang). Current: ${CMAKE_CXX_COMPILER}")
  endif()
endif()

# ---- Dependencies ----
find_package(ROOT 6.36 REQUIRED COMPONENTS
  Core RIO Hist Tree TreePlayer MathCore Graf Gpad Rint Thread 
)

# (Optional) Threads is safe to link on Linux; ROOT may bring it anyway
find_package(Threads REQUIRED)

# ---- Output dir for scripts ----
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/scripts)

# ---- Scripts scan (ignore macOS AppleDouble) ----
set(SCRIPTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/script)
file(GLOB SCRIPT_SOURCES CONFIGURE_DEPENDS
  ${SCRIPTS_DIR}/*.cxx
  ${SCRIPTS_DIR}/*.cpp
)
list(FILTER SCRIPT_SOURCES EXCLUDE REGEX "/\\._[^/]+\\.(cxx|cpp)$")

# ---- cxxopts via FetchContent ----
FetchContent_Declare(
  cxxopts
  GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
  GIT_TAG        v3.2.0
)
FetchContent_GetProperties(cxxopts)
if(NOT cxxopts_POPULATED)
  FetchContent_Populate(cxxopts)
  add_subdirectory(${cxxopts_SOURCE_DIR} ${cxxopts_BINARY_DIR})
endif()

# ---- spdlog via FetchContent ----
FetchContent_Declare(
  spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG        v1.14.1
)
FetchContent_MakeAvailable(spdlog)

# ---- binparse via FetchContent (your repo must export binparse::binparse) ----
FetchContent_Declare(
  binparse
  GIT_REPOSITORY https://github.com/RockPie/CRU_RDH_Decoder.git
  GIT_TAG        v0.5
)
FetchContent_MakeAvailable(binparse)

# ---- Build one exe per script ----
set(SCRIPT_TARGETS)
foreach(src IN LISTS SCRIPT_SOURCES)
  get_filename_component(exec_name ${src} NAME_WE)
  add_executable(${exec_name} ${src})

  target_include_directories(${exec_name} PRIVATE ${CMAKE_SOURCE_DIR}/include)

  target_link_libraries(${exec_name}
    PRIVATE
        cxxopts
        spdlog::spdlog
        binparse::binparse
        ROOT::Core ROOT::RIO ROOT::Hist ROOT::Tree ROOT::MathCore
        ROOT::Graf ROOT::Gpad
        Threads::Threads
        ROOT::RHTTP
  )

  target_compile_definitions(${exec_name} PRIVATE STAGE_NAME="${exec_name}")

  list(APPEND SCRIPT_TARGETS ${exec_name})
endforeach()

if(SCRIPT_TARGETS)
  add_custom_target(all_scripts DEPENDS ${SCRIPT_TARGETS})
endif()

# ---- Helpful diagnostics ----
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
if(APPLE)
  message(STATUS "macOS deployment target: ${CMAKE_OSX_DEPLOYMENT_TARGET}")
endif()