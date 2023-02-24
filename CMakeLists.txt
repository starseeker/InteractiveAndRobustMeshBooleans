cmake_minimum_required(VERSION 3.21)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set the project name
project(mesh_booleans)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if (NOT DEFINED BIN_DIR)
  set(BIN_DIR bin)
endif (NOT DEFINED BIN_DIR)
if (NOT DEFINED LIB_DIR)
  set(LIB_DIR lib)
endif (NOT DEFINED LIB_DIR)

set(TBB_TEST OFF CACHE BOOL " " FORCE)
set(TBB_EXAMPLES OFF CACHE BOOL " " FORCE)
add_subdirectory(arrangements/external/oneTBB)

set(cinolib_DIR ${PROJECT_SOURCE_DIR}/arrangements/external/Cinolib)
set(CINOLIB_USES_SHEWCHUK_PREDICATES ON)
set(CINOLIB_USES_INDIRECT_PREDICATES ON)

if (NOT DEFINED BUILD_EXAMPLES)
  set(BUILD_EXAMPLES ON)
endif (NOT DEFINED BUILD_EXAMPLES)
if (BUILD_EXAMPLES)
  set(CINOLIB_USES_OPENGL_GLFW_IMGUI ON)
endif (BUILD_EXAMPLES)

find_package(cinolib REQUIRED)

set(IRMB_INCLUDE_DIRS
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/code
  ${PROJECT_SOURCE_DIR}/arrangements/code
  ${PROJECT_SOURCE_DIR}/arrangements/external/abseil-cpp
  ${PROJECT_SOURCE_DIR}/arrangements/external/oneTBB
  )

set(IRMB_LIBS cinolib tbb)

# Build shared library
add_library(irmb SHARED irmb.cpp)
target_include_directories(irmb PUBLIC ${IRMB_INCLUDE_DIRS})
target_link_libraries(irmb ${IRMB_LIBS})
target_compile_definitions(irmb PUBLIC TBB_PARALLEL=1)
set_target_properties(irmb PROPERTIES VERSION 1.0.1 SOVERSION 1)
install(TARGETS irmb
  RUNTIME DESTINATION ${BIN_DIR}
  LIBRARY DESTINATION ${LIB_DIR}
  ARCHIVE DESTINATION ${LIB_DIR})
install(FILES irmb.h DESTINATION include/irmb)

if (BUILD_EXAMPLES)
  # add the executable
  add_executable(${PROJECT_NAME} main.cpp)
  target_include_directories(${PROJECT_NAME} PUBLIC ${IRMB_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} ${IRMB_LIBS})
  target_compile_definitions(${PROJECT_NAME} PUBLIC TBB_PARALLEL=1)

  add_executable(${PROJECT_NAME}_rotation main-rotation.cpp)
  target_include_directories(${PROJECT_NAME}_rotation PUBLIC ${IRMB_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME}_rotation ${IRMB_LIBS})
  target_compile_definitions(${PROJECT_NAME}_rotation PUBLIC TBB_PARALLEL=1)

  add_executable(${PROJECT_NAME}_arap main-arap.cpp)
  target_include_directories(${PROJECT_NAME}_arap PUBLIC ${IRMB_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME}_arap ${IRMB_LIBS})
  target_compile_definitions(${PROJECT_NAME}_arap PUBLIC TBB_PARALLEL=1)

  add_executable(${PROJECT_NAME}_stencil main-stencil.cpp)
  target_include_directories(${PROJECT_NAME}_stencil PUBLIC ${IRMB_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME}_stencil ${IRMB_LIBS})
  target_compile_definitions(${PROJECT_NAME}_stencil PUBLIC TBB_PARALLEL=1)

  add_executable(${PROJECT_NAME}_inputcheck main-inputcheck.cpp)
  target_include_directories(${PROJECT_NAME}_inputcheck PUBLIC ${IRMB_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME}_inputcheck ${IRMB_LIBS})
  target_compile_definitions(${PROJECT_NAME}_inputcheck PUBLIC TBB_PARALLEL=1)
endif (BUILD_EXAMPLES)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
