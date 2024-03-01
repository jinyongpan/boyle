add_header_only_library(math_curve2 curve2.hpp)
add_library(math::curve2 ALIAS math_curve2)

add_subdirectory(piecewise_curves)
