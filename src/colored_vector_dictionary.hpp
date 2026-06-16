#pragma once

#include "color.hpp"
#include "geometry.hpp"
#include "vector_key.hpp"

#include <cstddef>
#include <unordered_map>

namespace ot {

    struct colored_vector_key {
        vector_key vertex;
        color vertex_color;

        colored_vector_key() = default;
        colored_vector_key(double epsilon, const ot::vector3& value, ot::color color_);

        friend bool operator==(const colored_vector_key& lhs, const colored_vector_key& rhs) = default;
    };

    struct colored_vector_key_hash {
        std::size_t operator()(const colored_vector_key& key) const noexcept;
    };

    template <typename T>
    class colored_vector_dictionary {
    public:
        explicit colored_vector_dictionary(double eps)
            : epsilon_(2.0 * eps)
        {}

        bool contains(const ot::vector3& value, ot::color color) const
        {
            return values_.contains(make_key(value, color));
        }

        T& operator()(const ot::vector3& value, ot::color color)
        {
            return values_[make_key(value, color)];
        }

        const T& at(const ot::vector3& value, ot::color color) const
        {
            return values_.at(make_key(value, color));
        }

        T& at(const ot::vector3& value, ot::color color)
        {
            return values_.at(make_key(value, color));
        }

        void clear()
        {
            values_.clear();
        }

        void remove(const ot::vector3& value, ot::color color)
        {
            values_.erase(make_key(value, color));
        }

    private:
        colored_vector_key make_key(const ot::vector3& value, ot::color color) const
        {
            return colored_vector_key(epsilon_, value, color);
        }

        double epsilon_ = 0.0;
        std::unordered_map<colored_vector_key, T, colored_vector_key_hash> values_;
    };
}
