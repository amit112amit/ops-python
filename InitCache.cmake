# Use this to initialize the cache when building on Windows with MSYS Makefiles
# e.g. cmake -c InitCache.cmake -G "MSYS Makefiles" /path/to/ops-python

set(PYTHONLIBS_FOUND "ON" CACHE BOOL "pybind11 - have the Python libs been found")
set(PYTHON_PREFIX "C:/Program Files/Python37" CACHE PATH "pybind11 - path to the Python installation")
set(PYTHON_LIBRARIES "C:/Program Files/Python37/libs/python37.lib" CACHE FILEPATH "pybind11 - path to the python library")
set(PYTHON_INCLUDE_DIRS "C:/Program Files/Python37/include" CACHE PATH "pybind11 - path to where Python.h is found")
set(PYTHON_MODULE_EXTENSION ".pyd" CACHE STRING "pybind11 - lib extension, e.g. '.so' or '.pyd'")
set(PYTHON_MODULE_PREFIX "" CACHE STRING "pybind11 - lib name prefix: usually an empty string")
set(PYTHON_SITE_PACKAGES "C:/Program Files/Python37/Lib/site-packages" CACHE PATH "pybind11 - path to installation site-packages")
set(PYTHON_IS_DEBUG "OFF" CACHE BOOL "pybind11 - whether the Python interpreter is a debug build")