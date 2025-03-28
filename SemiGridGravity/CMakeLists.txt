file(
  DOWNLOAD
  https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.38.3/CPM.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
  EXPECTED_HASH SHA256=cc155ce02e7945e7b8967ddfaff0b050e958a723ef7aad3766d368940cb15494
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

CPMAddPackage(
    "gh:raysan5/raylib#5.5"
)

CPMAddPackage(
    NAME raylib-cpp
    GITHUB_REPOSITORY RobLoach/raylib-cpp
    VERSION 5.5.0
    OPTIONS 
        "BUILD_EXAMPLES OFF"
    DEPENDS raylib
)

add_executable (SemiGridGravity
	Types.hpp
	Particle.hpp
	Particle.cpp
    Simulation.hpp
	Simulation.cpp
	main.cpp
)

target_link_libraries(SemiGridGravity 
    raylib 
    raylib_cpp
)

set_property(TARGET SemiGridGravity PROPERTY CXX_STANDARD 23)
