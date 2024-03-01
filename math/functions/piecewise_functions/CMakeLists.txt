add_header_only_library(math_piecewise_linear_function1 piecewise_linear_function1.hpp)
target_link_libraries(math_piecewise_linear_function1
  INTERFACE
    Boost::serialization
    math::function1
    math::utils
)
add_library(math::piecewise_linear_function1 ALIAS math_piecewise_linear_function1)
add_doctest_test(math_piecewise_linear_function1 piecewise_linear_function1_test.cpp)

add_header_only_library(math_piecewise_cubic_function1 piecewise_cubic_function1.hpp)
target_link_libraries(math_piecewise_cubic_function1
  INTERFACE
    Boost::serialization
    math::function1
    math::cubic_interpolation
)
add_library(math::piecewise_cubic_function1 ALIAS math_piecewise_cubic_function1)
add_doctest_test(math_piecewise_cubic_function1 piecewise_cubic_function1_test.cpp)

add_header_only_library(math_piecewise_quintic_function1 piecewise_quintic_function1.hpp)
target_link_libraries(math_piecewise_quintic_function1
  INTERFACE
    Boost::serialization
    math::function1
    math::quintic_interpolation
)
add_library(math::piecewise_quintic_function1 ALIAS math_piecewise_quintic_function1)
add_doctest_test(math_piecewise_quintic_function1 piecewise_quintic_function1_test.cpp)
