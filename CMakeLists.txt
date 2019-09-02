cmake_minimum_required(VERSION 3.0.2)

project(OPS CXX Fortran)

cmake_policy(SET CMP0053 NEW)
cmake_policy(SET CMP0054 NEW)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Optimized compilation flags." FORCE)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(VTK REQUIRED)
vtk_module_config(VTK
	vtkCommonCore
	vtkCommonDataModel
	vtkFiltersCore
	vtkFiltersExtraction
	vtkFiltersGeometry
	vtkIOLegacy)

find_package(CGAL REQUIRED)

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

# A library to expose C-style function
add_library(opslib SHARED ExportConvexHull.cxx)
target_link_libraries(opslib PUBLIC CGAL::CGAL Eigen3::Eigen)

# Now create Python bindings
add_subdirectory(pybind11)
pybind11_add_module(ops OPSModelPython.cxx OPSModel.cxx LBFGSBWrapper.cxx ConvexHull.cxx)
target_include_directories(ops PUBLIC ${VTK_INCLUDE_DIRS})
target_link_libraries(ops PUBLIC CGAL::CGAL lbfgsb Eigen3::Eigen ${VTK_LIBRARIES})
target_compile_definitions(ops PUBLIC "_NO_PRINTING_" ${VTK_DEFINITIONS})
