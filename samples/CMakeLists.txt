add_library(gal_samples main.cpp)
target_link_libraries(gal_samples PRIVATE gal)

target_compile_options(gal_samples PRIVATE -save-temps=obj)