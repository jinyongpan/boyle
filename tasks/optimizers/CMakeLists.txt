add_library(route_line_builder STATIC route_line_builder.cpp)
target_sources(route_line_builder
  PUBLIC
    route_line_builder.h
  PRIVATE
    route_line_builder.cpp
)
target_link_libraries(route_line_builder
  INTERFACE
    task
    workspace
)

add_library(path_optimizer STATIC path_optimizer.cpp)
target_sources(path_optimizer
  PUBLIC
    path_optimizer.h
  PRIVATE
    path_optimizer.cpp
)
target_link_libraries(path_optimizer
  INTERFACE
    task
    workspace
)

add_library(speed_optimizer STATIC speed_optimizer.cpp)
target_sources(speed_optimizer
  PUBLIC
    speed_optimizer.h
  PRIVATE
    speed_optimizer.cpp
)
target_link_libraries(path_optimizer
  INTERFACE
    task
    workspace
)
