cmake_minimum_required(VERSION 3.0.2)

project(OPS CXX Fortran)

cmake_policy(SET CMP0053 NEW)
cmake_policy(SET CMP0054 NEW)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Optimized compilation flags." FORCE)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(CGAL REQUIRED)
include( ${CGAL_USE_FILE} )

find_package(Boost REQUIRED COMPONENTS "filesystem")
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include(CMakeAddFortranSubdirectory)
cmake_add_fortran_subdirectory( lbfgsb
		PROJECT lbfgsb
		ARCHIVE_DIR ../
		RUNTIME_DIR ../
		LIBRARIES lbfgsb
		LINK_LIBRARIES
		NO_EXTERNAL_INSTALL
		)

# Now create Python bindings
add_subdirectory(pybind11)
pybind11_add_module(ops OPSModelPython.cxx OPSModel.cxx LBFGSBWrapper.cxx)
target_link_libraries(ops PUBLIC ${VTK_LIBRARIES} ${CGAL_LIBRARIES} lbfgsb Eigen3::Eigen)
target_include_directories(ops PUBLIC ${VTK_INCLUDE_DIRS})
target_compile_definitions(ops PUBLIC "_NO_PRINTING_")
