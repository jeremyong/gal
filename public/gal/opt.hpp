#pragma once

#if defined(__clang) || defined(__GNUG__)
#    define GAL_FORCE_INLINE __attribute__((always_inline))
#elif defined(_MSC_VER)
#    define GAL_FORCE_INLINE __forceinline
#endif

#define GAL_NODISCARD [[nodiscard]]
