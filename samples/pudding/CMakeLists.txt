include(FetchContent)

FetchContent_Declare(
    sdl2
    GIT_REPOSITORY https://github.com/SDL-mirror/SDL.git
    GIT_TAG release-2.0.10
)

FetchContent_Declare(
    fmtlib
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG 6.0.0
)

FetchContent_Declare(
    vma
    GIT_REPOSITORY https://github.com/GPUOpen-LibrariesAndSDKs/VulkanMemoryAllocator.git
    GIT_TAG v2.2.0
)

FetchContent_MakeAvailable(sdl2 fmtlib vma)

# For now, link against installed SDK
find_package(Vulkan REQUIRED)

add_library(pudding
    src/cmd_stream.cpp
    src/device.cpp
    src/gc.cpp
    src/image.cpp
    src/image_view.cpp
    src/log.cpp
    src/pudding.cpp
    src/renderer.cpp
    src/swapchain.cpp
    src/window.cpp
    src/vk_mem_alloc.cpp
    src/uid.cpp
)

target_include_directories(pudding
    PRIVATE public/pudding
    PRIVATE ${vma_SOURCE_DIR}/src
    PUBLIC public
)

target_link_libraries(pudding
    PRIVATE gal
    PRIVATE Vulkan::Vulkan
    PRIVATE SDL2
    PRIVATE fmt::fmt
)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_compile_options(pudding PUBLIC -pipe)
endif()