add_library(kinetics_route_line_cubic_acc_model STATIC route_line_cubic_acc_model.cpp)
target_sources(kinetics_route_line_cubic_acc_model
  PUBLIC
    route_line_cubic_acc_model.h
  PRIVATE
    route_line_cubic_acc_model.cpp
)
target_link_libraries(kinetics_route_line_cubic_acc_model
  PUBLIC
    math::qp_problem
    math::osqp_solver
    kinetics::motion1
    kinetics::fence1
)
add_library(kinetics::route_line_cubic_acc_model ALIAS kinetics_route_line_cubic_acc_model)
add_doctest_test(kinetics_route_line_cubic_acc_model route_line_cubic_acc_model_test.cpp)

add_library(kinetics_route_line_quintic_acc_model STATIC route_line_quintic_acc_model.cpp)
target_sources(kinetics_route_line_quintic_acc_model
  PUBLIC
    route_line_quintic_acc_model.h
  PRIVATE
    route_line_quintic_acc_model.cpp
)
target_link_libraries(kinetics_route_line_quintic_acc_model
  PUBLIC
    math::qp_problem
    math::osqp_solver
    kinetics::motion1
    kinetics::fence1
)
add_library(kinetics::route_line_quintic_acc_model ALIAS kinetics_route_line_quintic_acc_model)
add_doctest_test(kinetics_route_line_quintic_acc_model route_line_quintic_acc_model_test.cpp)
