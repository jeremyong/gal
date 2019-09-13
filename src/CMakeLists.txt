add_library(gal INTERFACE)
target_include_directories(gal INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(gal INTERFACE fmt)