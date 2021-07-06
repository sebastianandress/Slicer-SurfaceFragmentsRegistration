cmake_minimum_required(VERSION 3.13.4)

project(Surface Fragments Registration)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://www.slicer.org/wiki/Documentation/Nightly/Extensions/SurfaceFragmentsRegistration")
set(EXTENSION_CATEGORY "Registration")
set(EXTENSION_CONTRIBUTORS "Sebastian Andress (LMU Munich)")
set(EXTENSION_DESCRIPTION " ")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/SurfaceFragmentsRegistration.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/SurfaceFragmentsRegistration/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(SurfaceFragmentsRegistration)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
