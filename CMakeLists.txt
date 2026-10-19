cmake_minimum_required(VERSION 3.25 FATAL_ERROR)

project(eic_event_generator)

# -----------------------------------------------------------------------------
# 1. Run sed to extract n_fdet and n_rdet from detector.in
# -----------------------------------------------------------------------------
set(DETECTOR_IN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/build/input/detector.in")

# Re-run CMake if detector.in changes
set_property(
    DIRECTORY 
    APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS 
    "${DETECTOR_IN_FILE}"
)

# Extract n_fdet value using sed
execute_process(
    COMMAND bash -c "sed -n 's/^n_fdet[[:space:]]*=[[:space:]]*\\([0-9]*\\)/\\1/p' '${DETECTOR_IN_FILE}'"
    OUTPUT_VARIABLE N_FDET
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Extract n_rdet value using sed
execute_process(
    COMMAND bash -c "sed -n 's/^n_rdet[[:space:]]*=[[:space:]]*\\([0-9]*\\)/\\1/p' '${DETECTOR_IN_FILE}'"
    OUTPUT_VARIABLE N_RDET
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "Sed extracted N_FDET = ${N_FDET}, N_RDET = ${N_RDET}")

# -----------------------------------------------------------------------------
# 2. Build Target Configuration
# -----------------------------------------------------------------------------
find_package(ROOT REQUIRED COMPONENTS Graf Gpad Hist)
include(${ROOT_USE_FILE})

file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB beam ${PROJECT_SOURCE_DIR}/beam/*.cc)
file(GLOB target ${PROJECT_SOURCE_DIR}/target/*.cc)
file(GLOB event ${PROJECT_SOURCE_DIR}/event/*.cc)
file(GLOB detector ${PROJECT_SOURCE_DIR}/detector/*.cc)

add_executable(eic_event_generator ${headers} ${sources} ${beam} ${target} ${event} ${detector} main.cpp)

# Pass N_FDET and N_RDET as separate macro definitions to C++
target_compile_definitions(eic_event_generator PRIVATE
    N_FDET=${N_FDET}
    N_RDET=${N_RDET}
)

install(TARGETS eic_event_generator RUNTIME DESTINATION bin)
target_link_libraries(eic_event_generator ROOT::Core ROOT::Graf ROOT::Gpad ROOT::Hist)
target_compile_features(eic_event_generator PRIVATE cxx_std_17)

include_directories(${PROJECT_SOURCE_DIR}/include)