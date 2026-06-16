#include "color.hpp"

#include <functional>

ot::color::color() = default;

ot::color::color(std::uint8_t red_, std::uint8_t green_, std::uint8_t blue_)
    : red(red_)
    , green(green_)
    , blue(blue_)
{}

ot::color::color(std::uint32_t rgb)
    : red(static_cast<std::uint8_t>((rgb >> 16) & 0xff))
    , green(static_cast<std::uint8_t>((rgb >> 8) & 0xff))
    , blue(static_cast<std::uint8_t>(rgb & 0xff))
{}

std::uint32_t ot::color::to_rgb() const
{
    std::uint32_t rgb = static_cast<std::uint32_t>(red);
    rgb = (rgb << 8) + static_cast<std::uint32_t>(green);
    rgb = (rgb << 8) + static_cast<std::uint32_t>(blue);
    return rgb;
}

std::size_t ot::color_hash::operator()(const ot::color& value) const noexcept
{
    std::size_t seed = 17;
    seed ^= std::hash<std::uint8_t>{}(value.red) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    seed ^= std::hash<std::uint8_t>{}(value.green) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    seed ^= std::hash<std::uint8_t>{}(value.blue) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    return seed;
}
