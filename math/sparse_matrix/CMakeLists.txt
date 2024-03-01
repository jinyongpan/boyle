add_header_only_library(math_sparse_matrix sparse_matrix.hpp)
add_library(math::sparse_matrix ALIAS math_sparse_matrix)

add_header_only_library(math_coo_matrix coo_matrix.hpp)
target_link_libraries(math_coo_matrix
  INTERFACE
    Boost::serialization
    math::sparse_matrix
)
add_library(math::coo_matrix ALIAS math_coo_matrix)
add_doctest_test(math_coo_matrix coo_matrix_test.cpp)

add_header_only_library(math_lil_matrix lil_matrix.hpp)
target_link_libraries(math_lil_matrix
  INTERFACE
    Boost::serialization
    math::sparse_matrix
)
add_library(math::lil_matrix ALIAS math_lil_matrix)
add_doctest_test(math_lil_matrix lil_matrix_test.cpp)

add_header_only_library(math_csc_matrix csc_matrix.hpp)
target_link_libraries(math_csc_matrix
  INTERFACE
    Boost::serialization
    math::sparse_matrix
)
add_library(math::csc_matrix ALIAS math_csc_matrix)
add_doctest_test(math_csc_matrix csc_matrix_test.cpp)

add_header_only_library(math_csr_matrix csr_matrix.hpp)
target_link_libraries(math_csr_matrix
  INTERFACE
    Boost::serialization
    math::sparse_matrix
)
add_library(math::csr_matrix ALIAS math_csr_matrix)
add_doctest_test(math_csr_matrix csr_matrix_test.cpp)
