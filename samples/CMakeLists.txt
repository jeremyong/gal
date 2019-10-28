add_subdirectory(pudding)

add_executable(gal_samples main.cpp)
target_link_libraries(gal_samples PRIVATE gal pudding)

target_compile_definitions(gal_samples PRIVATE GAL_DEBUG)
# target_compile_options(gal_samples PRIVATE -save-temps=obj -O3 -DNDEBUG)
# target_compile_options(gal_samples PRIVATE -ftime-report -ftime-report-details -Q)
# target_compile_options(gal_samples PRIVATE -ftime-trace)