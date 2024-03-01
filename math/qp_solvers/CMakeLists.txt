add_header_only_library(math_qp_problem qp_problem.hpp)
target_link_libraries(math_qp_problem
  INTERFACE
    Boost::serialization
    math::coo_matrix
    math::lil_matrix
    math::utils
)
add_library(math::qp_problem ALIAS math_qp_problem)
add_doctest_test(math_qp_problem qp_problem_test.cpp)

add_header_only_library(math_osqp_solver osqp_solver.hpp)
target_link_libraries(math_osqp_solver
  INTERFACE
    Boost::serialization
    osqpstatic
    math::csc_matrix
    math::qp_problem
)
add_library(math::osqp_solver ALIAS math_osqp_solver)
add_doctest_test(math_osqp_solver osqp_solver_test.cpp)
