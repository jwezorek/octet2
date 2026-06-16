#pragma once

#include "geometry.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ot {

struct vector_key {
    std::int64_t x = 0;
    std::int64_t y = 0;
    std::int64_t z = 0;

    vector_key() = default;

    vector_key(double epsilon, const ot::vector3& value)
        : x(float_to_key_value(value.x(), epsilon))
        , y(float_to_key_value(value.y(), epsilon))
        , z(float_to_key_value(value.z(), epsilon))
    {
    }

    [[nodiscard]] static std::int64_t float_to_key_value(double value, double epsilon)
    {
        return static_cast<std::int64_t>(std::llround(value / epsilon));
    }

    friend bool operator==(const vector_key& lhs, const vector_key& rhs) = default;
    friend auto operator<=>(const vector_key& lhs, const vector_key& rhs) = default;
};

struct vector_key_hash {
    [[nodiscard]] std::size_t operator()(const vector_key& key) const noexcept
    {
        std::size_t seed = 17;
        combine(seed, std::hash<std::int64_t>{}(key.x));
        combine(seed, std::hash<std::int64_t>{}(key.y));
        combine(seed, std::hash<std::int64_t>{}(key.z));
        return seed;
    }

private:
    static void combine(std::size_t& seed, std::size_t value) noexcept
    {
        seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    }
};

template <typename T>
class vector_dictionary {
public:
    explicit vector_dictionary(double eps)
        : epsilon_(2.0 * eps)
    {
    }

    [[nodiscard]] bool contains(const ot::vector3& value) const
    {
        return values_.contains(make_key(value));
    }

    T& operator[](const ot::vector3& value)
    {
        return values_[make_key(value)];
    }

    const T& at(const ot::vector3& value) const
    {
        return values_.at(make_key(value));
    }

    T& at(const ot::vector3& value)
    {
        return values_.at(make_key(value));
    }

    void clear()
    {
        values_.clear();
    }

    void remove(const ot::vector3& value)
    {
        values_.erase(make_key(value));
    }

private:
    [[nodiscard]] vector_key make_key(const ot::vector3& value) const
    {
        return vector_key(epsilon_, value);
    }

    double epsilon_ = 0.0;
    std::unordered_map<vector_key, T, vector_key_hash> values_;
};

class vector_set {
public:
    explicit vector_set(double epsilon)
        : values_(epsilon)
    {
    }

    void clear()
    {
        values_.clear();
    }

    void add(const ot::vector3& value)
    {
        values_[value] = true;
    }

    void add_range(const std::vector<ot::vector3>& values)
    {
        for (const ot::vector3& value : values) {
            add(value);
        }
    }

    [[nodiscard]] bool contains(const ot::vector3& value) const
    {
        return values_.contains(value);
    }

private:
    vector_dictionary<bool> values_;
};

struct vec_list_key {
    std::vector<vector_key> keys;

    vec_list_key() = default;

    vec_list_key(double epsilon, const std::vector<ot::vector3>& values, bool order_is_significant)
    {
        keys.reserve(values.size());
        for (const ot::vector3& value : values) {
            keys.emplace_back(epsilon, value);
        }

        if (!order_is_significant) {
            std::ranges::sort(keys);
        }
    }

    friend bool operator==(const vec_list_key& lhs, const vec_list_key& rhs) = default;
};

struct vec_list_key_hash {
    [[nodiscard]] std::size_t operator()(const vec_list_key& key) const noexcept
    {
        std::size_t seed = 17;
        vector_key_hash hasher;
        for (const vector_key& value : key.keys) {
            seed ^= hasher(value) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

template <typename T>
class vec_list_dictionary {
public:
    vec_list_dictionary(double eps, bool order_is_significant)
        : epsilon_(2.0 * eps)
        , order_is_significant_(order_is_significant)
    {
    }

    T& operator[](const std::vector<ot::vector3>& values)
    {
        return values_[make_key(values)];
    }

    void insert_or_assign(const std::vector<ot::vector3>& values, const T& value)
    {
        values_.insert_or_assign(make_key(values), value);
    }

    void insert_or_assign(const std::vector<ot::vector3>& values, T&& value)
    {
        values_.insert_or_assign(make_key(values), std::move(value));
    }

    const T& at(const std::vector<ot::vector3>& values) const
    {
        return values_.at(make_key(values));
    }

    [[nodiscard]] bool contains(const std::vector<ot::vector3>& values) const
    {
        return values_.contains(make_key(values));
    }

    void clear()
    {
        values_.clear();
    }

    void remove(const std::vector<ot::vector3>& values)
    {
        values_.erase(make_key(values));
    }

    [[nodiscard]] std::vector<T> values() const
    {
        std::vector<T> output;
        output.reserve(values_.size());
        for (const auto& [key, value] : values_) {
            output.push_back(value);
        }
        return output;
    }

    [[nodiscard]] std::size_t size() const
    {
        return values_.size();
    }

private:
    [[nodiscard]] vec_list_key make_key(const std::vector<ot::vector3>& values) const
    {
        return vec_list_key(epsilon_, values, order_is_significant_);
    }

    double epsilon_ = 0.0;
    bool order_is_significant_ = true;
    std::unordered_map<vec_list_key, T, vec_list_key_hash> values_;
};

class vec_list_set {
public:
    vec_list_set(double epsilon, bool order_is_significant)
        : values_(epsilon, order_is_significant)
    {
    }

    void clear()
    {
        values_.clear();
    }

    void add(const std::vector<ot::vector3>& values)
    {
        values_[values] = true;
    }

    void add_range(const std::vector<std::vector<ot::vector3>>& values)
    {
        for (const auto& value : values) {
            add(value);
        }
    }

    [[nodiscard]] bool contains(const std::vector<ot::vector3>& values) const
    {
        return values_.contains(values);
    }

    [[nodiscard]] std::size_t size() const
    {
        return values_.size();
    }

private:
    vec_list_dictionary<bool> values_;
};

} // namespace octet
