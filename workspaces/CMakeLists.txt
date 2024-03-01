add_header_only_library(workspace workspace.hpp)

add_header_only_library(common_workspace common_workspace.hpp)
target_link_libraries(common_workspace
  INTERFACE
    workspace
    math::vec2
)
add_doctest_test(common_workspace common_workspace_test.cpp)

add_header_only_library(config_workspace config_workspace.hpp)

add_header_only_library(lateral_workspace lateral_workspace.hpp)

add_header_only_library(longitudinal_workspace longitudinal_workspace.hpp)

add_header_only_library(prediction_workspace_test predition_workspace.hpp)

add_header_only_library(signal_workspace signal_workspace.hpp)
