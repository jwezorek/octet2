#pragma once

#include "cell.hpp"
#include "octet.hpp"

#include <cstddef>
#include <optional>
#include <unordered_set>
#include <vector>

namespace ot {

    class cell_set {
    public:
        cell_set();

        explicit cell_set(std::optional<cell> value);
        explicit cell_set(const cell& value);
        explicit cell_set(const std::vector<cell>& values);

        std::vector<cell> cells() const;

        std::vector<cell> octahedra() const;
        std::vector<cell> tetrahedra() const;
        std::vector<cell> pyramids() const;
        std::vector<cell> octahedra_and_pyramids() const;

        bool has_cell(const cell& value) const;

        void add_cell(const cell& value);
        void remove_cell(const cell& value);

        void add_cells(const std::vector<cell>& values);
        void add_cells(const cell_set& values);

        cell_set intersect(const cell_set& other) const;
        cell_set intersect(const std::vector<cell>& values) const;

        cell_set unite(const cell_set& other) const;
        cell_set unite(const std::vector<cell>& values) const;

        cell_set subtract(const cell_set& other) const;

        cell_set translate(int column_delta, int row_delta, int layer_delta) const;
        cell_set apply(const ot::matrix4& matrix) const;

        int pyramid_weight() const;
        int weight() const;

        std::size_t count() const;
        bool any() const;

        std::optional<cell> contents_of_octahedral_cell(const cell& index) const;

        bool can_add(const cell_set& values) const;

    private:
        std::vector<cell> select_by_type(cell_type type) const;

        std::unordered_set<cell, cell_hash> cells_;
    };

} // namespace ot