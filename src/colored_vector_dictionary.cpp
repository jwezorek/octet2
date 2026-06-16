#include "colored_vector_dictionary.hpp"

ot::colored_vector_key::colored_vector_key(double epsilon, const ot::vector3& value, ot::color color_)
    : vertex(epsilon, value)
    , vertex_color(color_)
{}

std::size_t ot::colored_vector_key_hash::operator()(const ot::colored_vector_key& key) const noexcept
{
    std::size_t seed = 17;
    seed ^= ot::vector_key_hash{}(key.vertex) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    seed ^= ot::color_hash{}(key.vertex_color) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    return seed;
}
