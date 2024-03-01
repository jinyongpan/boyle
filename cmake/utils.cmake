function(add_header_only_library name header_name)
  add_library(${name} INTERFACE)
  target_compile_definitions(${name} INTERFACE LIBRARY_HEADER_ONLY)
  target_sources(${name}
    INTERFACE
      ${header_name}
  )
endfunction()

function(add_doctest_test name source_name)
  set(test_name ${name}_test)
  add_executable(${test_name} ${source_name})
  target_link_libraries(${test_name}
    PRIVATE
      cxxopts::cxxopts
      doctest::doctest
      Matplot++::matplot
      spdlog::spdlog
      ${name}
  )
  add_test(
    NAME ${test_name}
    COMMAND ${test_name}
  )
endfunction()
