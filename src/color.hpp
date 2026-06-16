#pragma once

#include <cstddef>
#include <cstdint>

namespace ot {

    struct color {
        std::uint8_t red = 0;
        std::uint8_t green = 0;
        std::uint8_t blue = 0;

        color();
        color(std::uint8_t red_, std::uint8_t green_, std::uint8_t blue_);
        explicit color(std::uint32_t rgb);

        std::uint32_t to_rgb() const;

        friend bool operator==(const color& lhs, const color& rhs) = default;
    };

    struct color_hash {
        std::size_t operator()(const color& value) const noexcept;
    };
}
