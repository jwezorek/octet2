#pragma once

#include "cell.hpp"

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ot {

    template <typename T>
    class colored_cell_set {
    public:
        colored_cell_set() = default;

        explicit colored_cell_set(std::optional<std::pair<cell, T>> value)
        {
            if (value) {
                add_cell(value->first, value->second);
            }
        }

        explicit colored_cell_set(const std::vector<std::pair<cell, T>>& values)
        {
            add_cells(values);
        }

        std::vector<std::pair<cell, T>> colored_cells() const
        {
            std::vector<std::pair<cell, T>> output;
            output.reserve(cells_.size());
            for (const auto& [value, color] : cells_) {
                output.emplace_back(value, color);
            }
            return output;
        }

        std::vector<cell> cells() const
        {
            std::vector<cell> output;
            output.reserve(cells_.size());
            for (const auto& [value, unused] : cells_) {
                output.push_back(value);
            }
            return output;
        }

        std::vector<std::pair<cell, T>> octahedra() const
        {
            return select_by_type(cell_type::octahedron);
        }

        std::vector<std::pair<cell, T>> tetrahedra() const
        {
            return select_by_type(cell_type::tetrahedron);
        }

        std::vector<std::pair<cell, T>> pyramids() const
        {
            return select_by_type(cell_type::pyramid);
        }

        bool has_cell(const cell& value) const
        {
            return cells_.contains(value);
        }

        std::optional<T> get_cell_color(const cell& value) const
        {
            auto iter = cells_.find(value);
            if (iter == cells_.end()) {
                return std::nullopt;
            }
            return iter->second;
        }

        void add_cell(const cell& value, const T& color)
        {
            cells_.insert_or_assign(value, color);
        }

        void remove_cell(const cell& value)
        {
            cells_.erase(value);
        }

        void add_cells(const std::vector<std::pair<cell, T>>& values)
        {
            for (const auto& [value, color] : values) {
                add_cell(value, color);
            }
        }

        void add_cells(const std::vector<cell>& values, const T& color)
        {
            for (const cell& value : values) {
                add_cell(value, color);
            }
        }

        void add_cells(const colored_cell_set<T>& values)
        {
            add_cells(values.colored_cells());
        }

        colored_cell_set<T> intersect(const colored_cell_set<T>& other) const
        {
            colored_cell_set<T> output;
            for (const auto& [value, color] : cells_) {
                auto iter = other.cells_.find(value);
                if (iter != other.cells_.end() && iter->second == color) {
                    output.add_cell(value, color);
                }
            }
            return output;
        }

        colored_cell_set<T> unite(const colored_cell_set<T>& other) const
        {
            colored_cell_set<T> output(*this);
            output.add_cells(other);
            return output;
        }

        colored_cell_set<T> subtract(const colored_cell_set<T>& other) const
        {
            colored_cell_set<T> output;
            for (const auto& [value, color] : cells_) {
                auto iter = other.cells_.find(value);
                if (iter == other.cells_.end() || iter->second != color) {
                    output.add_cell(value, color);
                }
            }
            return output;
        }

        colored_cell_set<T> translate(int column_delta, int row_delta, int layer_delta) const
        {
            colored_cell_set<T> output;
            for (const auto& [value, color] : cells_) {
                output.add_cell(value.translate(column_delta, row_delta, layer_delta), color);
            }
            return output;
        }

        colored_cell_set<T> apply(const ot::matrix4& matrix) const
        {
            colored_cell_set<T> output;
            for (const auto& [value, color] : cells_) {
                output.add_cell(value.apply(matrix), color);
            }
            return output;
        }

        std::size_t count() const
        {
            return cells_.size();
        }

        bool any() const
        {
            return !cells_.empty();
        }

        std::optional<cell> contents_of_octahedral_cell(const cell& index) const
        {
            if (index.type != cell_type::octahedron) {
                throw std::invalid_argument("Expected an octahedron");
            }

            if (cells_.contains(index)) {
                return index;
            }

            for (octa_vert vert : ot::octa_verts()) {
                const cell pyramid = index.to_pyramid(vert);
                if (cells_.contains(pyramid)) {
                    return pyramid;
                }
            }

            return std::nullopt;
        }

    private:
        std::vector<std::pair<cell, T>> select_by_type(cell_type type) const
        {
            std::vector<std::pair<cell, T>> output;
            for (const auto& [value, color] : cells_) {
                if (value.type == type) {
                    output.emplace_back(value, color);
                }
            }
            return output;
        }

        std::unordered_map<cell, T, cell_hash> cells_;
    };
}
