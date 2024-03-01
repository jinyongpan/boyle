add_header_only_library(math_function1 function1.hpp)
add_library(math::function1 ALIAS math_function1)

add_subdirectory(piecewise_functions)
