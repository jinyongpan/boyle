add_header_only_library(common_array_view array_view.hpp)
target_link_libraries(common_array_view
  INTERFACE
    Eigen3::Eigen
)
add_library(common::array_view ALIAS common_array_view)
add_doctest_test(common_array_view array_view_test.cpp)

add_header_only_library(common_fsm fsm.hpp)
add_library(common::fsm ALIAS common_fsm)
add_doctest_test(common_fsm fsm_test.cpp)
