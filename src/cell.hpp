#pragma once
#include "octet.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace ot {

    class cell;
    class cell_face;

    cell operator+(const cell& lhs, const cell& rhs);
    bool operator==(const cell& lhs, const cell& rhs);
    bool operator==(const cell_face& lhs, const cell_face& rhs);

    class cell {
    public:
        int column = 0;
        int row = 0;
        int layer = 0;
        cell_type type = cell_type::octahedron;

        cell();
        cell(int column_, int row_, int layer_);
        cell(tetrahedron_orientation orientation, int column_, int row_, int layer_);
        cell(octa_vert pyramid_tip, int column_, int row_, int layer_);

        octa_vert pyramid_tip() const;
        tetrahedron_orientation tetra_orientation() const;
        cell to_octahedron() const;
        cell to_tetrahedron(tetrahedron_orientation orientation) const;
        cell to_pyramid(octa_vert vert) const;

        std::vector<cell> neighbors() const;
        std::vector<cell> neighbors_with_pyramids() const;
        std::vector<cell> edge_neighbors() const;
        std::vector<cell> edge_neighbors_with_pyramids() const;
        cell translate(int column_delta, int row_delta, int layer_delta) const;
        cell apply(const ot::matrix4& matrix) const;
        std::vector<cell_face> faces() const;
        std::vector<cell> vertex_indices() const;
        ot::vector3 center() const;

        friend cell operator+(const cell& lhs, const cell& rhs);
        friend bool operator==(const cell& lhs, const cell& rhs);

    private:
        cell apply_matrix_to_tetrahedron(const ot::matrix4& matrix) const;
        cell apply_matrix_to_pyramid(const ot::matrix4& matrix) const;

        std::optional<tetrahedron_orientation> tetra_orientation_;
        std::optional<octa_vert> pyramid_tip_;
    };

    struct cell_hash {
        std::size_t operator()(const cell& value) const noexcept;

    private:
        static void combine(std::size_t& seed, std::size_t value) noexcept;
    };

    class cell_face {
    public:
        class edge;

        explicit cell_face(cell value, octa_face face);
        explicit cell_face(cell value, tetra_face_or_vert face);
        explicit cell_face(cell value);

        const cell& owner_cell() const;
        bool is_square_face() const;
        bool is_from_tetrahedron() const;

        std::vector<cell> vertex_indices() const;
        std::vector<edge> edges() const;
        std::vector<ot::vector3> vertices() const;

        friend bool operator==(const cell_face& lhs, const cell_face& rhs);

    private:
        cell cell_;
        std::optional<octa_face> octa_face_;
        std::optional<tetra_face_or_vert> tetra_face_;
        bool is_square_face_ = false;
    };

    class cell_face::edge {
    public:
        cell_face face;
        int from = 0;
        int to = 0;

        edge(cell_face face_, int from_, int to_);

        std::vector<ot::vector3> to_vertices() const;
        ot::vector3 from_vertex() const;
        ot::vector3 to_vertex() const;
    };

    struct cell_face_hash {
        std::size_t operator()(const cell_face& face) const noexcept;

    private:
        static void combine(std::size_t& seed, std::size_t value) noexcept;
    };
}
