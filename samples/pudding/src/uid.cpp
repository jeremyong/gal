#include <detail/uid.hpp>

#include <functional>
#include <vector>

using namespace pd::detail;

// NOTE, the zero id is reserved
static unique_id::id_t next = 1u;
static std::vector<unique_id::id_t> free_ids;

unique_id::unique_id()
{
    if (free_ids.empty())
    {
        id = next++;
    }
    else
    {
        id = free_ids.back();
        free_ids.pop_back();
    }
    hash = std::hash<id_t>{}(id);
}

unique_id::unique_id(unique_id&& other)
    : id{other.id}
    , hash{other.hash}
{
    other.id = 0;
}

unique_id& unique_id::operator=(unique_id&& other)
{
    std::swap(id, other.id);
    std::swap(hash, other.hash);
    return *this;
}

void unique_id::reset()
{
    unique_id other;
    release();
    *this = std::move(other);
}

void unique_id::release()
{
    if (id != 0)
    {
        free_ids.emplace_back(id);
        id = 0;
    }
}

unique_id::~unique_id()
{
    release();
}