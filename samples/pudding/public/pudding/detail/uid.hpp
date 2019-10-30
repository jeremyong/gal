#pragma once

#include <cstddef>
#include <cstdint>

namespace pd
{
    namespace detail
    {
        // Exhausting all available handles will result in an exception.
        struct unique_id final
        {
            using id_t = uint32_t;
            using hash_t = size_t;

            unique_id();
            unique_id(unique_id&& other);
            unique_id& operator=(unique_id&& other);
            unique_id(unique_id const&) = delete;
            unique_id& operator=(unique_id const&) = delete;

            // Acquires a new id and hash in place, releasing the old identifier
            void reset();

            // Returns a previously leased identifier to the pool. This should be done as soon as the referenced
            // resource is reclaimed.
            void release();

            // If the identifier was not already released, it will get released at this point
            ~unique_id();

            id_t id;
            hash_t hash;
        };
    }
}