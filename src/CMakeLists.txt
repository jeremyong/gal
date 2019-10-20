add_library(gal
    INTERFACE
)
target_include_directories(gal
    INTERFACE ${PROJECT_SOURCE_DIR}/public
    )
if (GAL_SAMPLES_ENABLED)
  target_link_libraries(gal INTERFACE fmt)
endif()
target_compile_features(gal INTERFACE cxx_std_17)
