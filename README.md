I have another repository `oriented-particles` which is in C++. Here I generate Python bindings for the a subset of the C++ code. The goal is to learn making Python bindings. It is also helpful to have the core libraries in fast C++ and to write drivers in Python using `h5py`, `pandas` etc. for managing input/output of the program.

## Building on Windows 10 with MSYS2:
1. Install MSYS2 and then using the MSYS2 shell install the following packages:
    1. mingw-w64-x86_64-toolchain
    2. mingw-w64-x86_64-vtk
    3. mingw-w64-x86_64-cgal
    4. mingw-w64-x86_64-eigen3
2. If you have installed MSYS2 as per defaults then your `MSYS2_ROOT=C:\msys64`. You need to fix the hardcoded file paths in the following files:
    1. `${MSYS2_ROOT}/mingw64/lib/cmake/Qt5Gui/Qt5GuiConfixExtras.cmake`: Look for line starting with `_qt5gui_find_extra_libs`.
    2. `${MSYS2_ROOT}/mingw64/lib/cmake/CGAL/CGALExports.cmake`: Look for line starting with `set_target_properties`.
3. You will also need to use the CMake initial cache file provided in the root directory of this project `InitCacheMSYS2.cmake` by passing it to CMake e.g. `cmake -c InitCacheMSYS2.cmake /path/to/ops-python`.

## Building on Windows 10 with WSL (Ubuntu 18.04):
1. You will have to build CGAL from source as the release in the repository is 4.11 and we need 4.14.
2. You will have to pass a flag to CMake to be able to link the static Fortran solver library built in this project with the final dynamic Python module library. You should invoke CMake as `cmake -DCMAKE_Fortran_FLAGS="-fPIC" /path/to/ops-python`.

## Building on Linux:
1. It should be relatively straightforward if the latest version of CGAL is installed. 
    
