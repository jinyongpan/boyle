add_header_only_library(task task.hpp)
target_link_libraries(task
  INTERFACE
    workspace
)

add_subdirectory(optimizers)
