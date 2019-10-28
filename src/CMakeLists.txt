# This is a placeholder folder for if/when GAL no longer becomes an interface target
add_library(gal
    INTERFACE
)
target_include_directories(gal
    INTERFACE ${PROJECT_SOURCE_DIR}/public
    )
target_compile_features(gal INTERFACE cxx_std_17)
