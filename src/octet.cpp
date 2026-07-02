#include "octet.hpp"
#include "cell_set.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

    constexpr double half_root_2 = 0.707106781186547524400844362104849039;

    struct enum_pair_hash {
        template <typename A, typename B>
        std::size_t operator()(const std::pair<A, B>& value) const noexcept
        {
            std::size_t seed = 17;
            seed ^= std::hash<int>{}(ot::enum_index(value.first)) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
            seed ^= std::hash<int>{}(ot::enum_index(value.second)) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
            return seed;
        }
    };

    struct slot_order_hash {
        std::size_t operator()(const std::pair<ot::slot, int>& value) const noexcept
        {
            std::size_t seed = 17;
            seed ^= std::hash<int>{}(ot::enum_index(value.first)) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
            seed ^= std::hash<int>{}(value.second) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
            return seed;
        }
    };

    struct bool_order_hash {
        std::size_t operator()(const std::pair<bool, int>& value) const noexcept
        {
            std::size_t seed = 17;
            seed ^= std::hash<bool>{}(value.first) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
            seed ^= std::hash<int>{}(value.second) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
            return seed;
        }
    };

    using octa_face_map = std::unordered_map<ot::octa_face, ot::cell>;
    using tetra_neighbor_key = std::pair<ot::tetrahedron_orientation, ot::tetra_face_or_vert>;
    using tetra_neighbor_map = std::unordered_map<tetra_neighbor_key, ot::cell, enum_pair_hash>;
    using octa_edge_vert_map = std::unordered_map<ot::octa_edge, std::pair<ot::octa_vert, ot::octa_vert>>;
    using opposite_vert_map = std::unordered_map<ot::octa_vert, ot::octa_vert>;

    const std::unordered_map<ot::octa_vert, std::vector<ot::octa_face>>& octa_faces_around_vert_table()
    {
        static const std::unordered_map<ot::octa_vert, std::vector<ot::octa_face>> table = {
            { ot::octa_vert::top, { ot::octa_face::north_top, ot::octa_face::east_top, ot::octa_face::south_top, ot::octa_face::west_top } },
            { ot::octa_vert::north_east, { ot::octa_face::north_top, ot::octa_face::north_bottom, ot::octa_face::east_bottom, ot::octa_face::east_top } },
            { ot::octa_vert::south_east, { ot::octa_face::south_bottom, ot::octa_face::south_top, ot::octa_face::east_top, ot::octa_face::east_bottom } },
            { ot::octa_vert::south_west, { ot::octa_face::west_bottom, ot::octa_face::west_top, ot::octa_face::south_top, ot::octa_face::south_bottom } },
            { ot::octa_vert::north_west, { ot::octa_face::north_bottom, ot::octa_face::north_top, ot::octa_face::west_bottom, ot::octa_face::west_top } },
            { ot::octa_vert::bottom, { ot::octa_face::north_bottom, ot::octa_face::west_bottom, ot::octa_face::south_bottom, ot::octa_face::east_bottom } }
        };
        return table;
    }

    const std::unordered_map<ot::octa_vert, std::vector<ot::octa_vert>>& octa_base_of_pyramid_table()
    {
        static const std::unordered_map<ot::octa_vert, std::vector<ot::octa_vert>> table = {
            { ot::octa_vert::top, { ot::octa_vert::north_east, ot::octa_vert::north_west, ot::octa_vert::south_west, ot::octa_vert::south_east } },
            { ot::octa_vert::bottom, { ot::octa_vert::south_east, ot::octa_vert::south_west, ot::octa_vert::north_west, ot::octa_vert::north_east } },
            { ot::octa_vert::north_east, { ot::octa_vert::south_east, ot::octa_vert::bottom, ot::octa_vert::north_west, ot::octa_vert::top } },
            { ot::octa_vert::south_west, { ot::octa_vert::top, ot::octa_vert::north_west, ot::octa_vert::bottom, ot::octa_vert::south_east } },
            { ot::octa_vert::north_west, { ot::octa_vert::north_east, ot::octa_vert::bottom, ot::octa_vert::south_west, ot::octa_vert::top } },
            { ot::octa_vert::south_east, { ot::octa_vert::top, ot::octa_vert::south_west, ot::octa_vert::bottom, ot::octa_vert::north_east } }
        };
        return table;
    }

    const std::unordered_map<ot::octa_edge, ot::cell>& octa_edge_to_neighbor_table()
    {
        static const std::unordered_map<ot::octa_edge, ot::cell> table = {
            { ot::octa_edge::north_west_top, ot::cell(-1, 0, 1) },
            { ot::octa_edge::north_east_top, ot::cell(0, 0, 1) },
            { ot::octa_edge::south_east_top, ot::cell(0, -1, 1) },
            { ot::octa_edge::south_west_top, ot::cell(-1, -1, 1) },
            { ot::octa_edge::south_west_bottom, ot::cell(0, 0, -1) },
            { ot::octa_edge::south_east_bottom, ot::cell(1, 0, -1) },
            { ot::octa_edge::north_west_bottom, ot::cell(0, 1, -1) },
            { ot::octa_edge::north_east_bottom, ot::cell(1, 1, -1) },
            { ot::octa_edge::south, ot::cell(0, -1, 0) },
            { ot::octa_edge::east, ot::cell(1, 0, 0) },
            { ot::octa_edge::north, ot::cell(0, 1, 0) },
            { ot::octa_edge::west, ot::cell(-1, 0, 0) }
        };
        return table;
    }

    const octa_edge_vert_map& octa_edge_to_verts_table()
    {
        static const octa_edge_vert_map table = {
            { ot::octa_edge::north, { ot::octa_vert::north_east, ot::octa_vert::north_west } },
            { ot::octa_edge::west, { ot::octa_vert::north_west, ot::octa_vert::south_west } },
            { ot::octa_edge::south, { ot::octa_vert::south_west, ot::octa_vert::south_east } },
            { ot::octa_edge::east, { ot::octa_vert::south_east, ot::octa_vert::north_east } },
            { ot::octa_edge::north_east_top, { ot::octa_vert::top, ot::octa_vert::north_east } },
            { ot::octa_edge::south_east_top, { ot::octa_vert::top, ot::octa_vert::south_east } },
            { ot::octa_edge::south_west_top, { ot::octa_vert::top, ot::octa_vert::south_west } },
            { ot::octa_edge::north_west_top, { ot::octa_vert::top, ot::octa_vert::north_west } },
            { ot::octa_edge::north_east_bottom, { ot::octa_vert::bottom, ot::octa_vert::north_east } },
            { ot::octa_edge::north_west_bottom, { ot::octa_vert::bottom, ot::octa_vert::north_west } },
            { ot::octa_edge::south_west_bottom, { ot::octa_vert::bottom, ot::octa_vert::south_west } },
            { ot::octa_edge::south_east_bottom, { ot::octa_vert::bottom, ot::octa_vert::south_east } }
        };
        return table;
    }

    const octa_face_map& neighbor_offset_table()
    {
        static const octa_face_map table = {
            { ot::octa_face::west_top, ot::cell(ot::tetrahedron_orientation::horizontal, 0, 0, 0) },
            { ot::octa_face::north_top, ot::cell(ot::tetrahedron_orientation::vertical, 0, 1, 0) },
            { ot::octa_face::east_top, ot::cell(ot::tetrahedron_orientation::horizontal, 1, 0, 0) },
            { ot::octa_face::south_top, ot::cell(ot::tetrahedron_orientation::vertical, 0, 0, 0) },
            { ot::octa_face::west_bottom, ot::cell(ot::tetrahedron_orientation::vertical, 0, 1, -1) },
            { ot::octa_face::north_bottom, ot::cell(ot::tetrahedron_orientation::horizontal, 1, 1, -1) },
            { ot::octa_face::east_bottom, ot::cell(ot::tetrahedron_orientation::vertical, 1, 1, -1) },
            { ot::octa_face::south_bottom, ot::cell(ot::tetrahedron_orientation::horizontal, 1, 0, -1) }
        };
        return table;
    }

    const tetra_neighbor_map& tetra_neighbor_offset_table()
    {
        static const tetra_neighbor_map table = {
            { { ot::tetrahedron_orientation::horizontal, ot::tetra_face_or_vert::north }, ot::cell(-1, 0, 1) },
            { { ot::tetrahedron_orientation::horizontal, ot::tetra_face_or_vert::east }, ot::cell(0, 0, 0) },
            { { ot::tetrahedron_orientation::horizontal, ot::tetra_face_or_vert::south }, ot::cell(-1, -1, 1) },
            { { ot::tetrahedron_orientation::horizontal, ot::tetra_face_or_vert::west }, ot::cell(-1, 0, 0) },
            { { ot::tetrahedron_orientation::vertical, ot::tetra_face_or_vert::north }, ot::cell(0, 0, 0) },
            { { ot::tetrahedron_orientation::vertical, ot::tetra_face_or_vert::east }, ot::cell(0, -1, 1) },
            { { ot::tetrahedron_orientation::vertical, ot::tetra_face_or_vert::south }, ot::cell(0, -1, 0) },
            { { ot::tetrahedron_orientation::vertical, ot::tetra_face_or_vert::west }, ot::cell(-1, -1, 1) }
        };
        return table;
    }

    const opposite_vert_map& opposite_octa_vert_table()
    {
        static const opposite_vert_map table = {
            { ot::octa_vert::top, ot::octa_vert::bottom },
            { ot::octa_vert::bottom, ot::octa_vert::top },
            { ot::octa_vert::north_east, ot::octa_vert::south_west },
            { ot::octa_vert::south_west, ot::octa_vert::north_east },
            { ot::octa_vert::north_west, ot::octa_vert::south_east },
            { ot::octa_vert::south_east, ot::octa_vert::north_west }
        };
        return table;
    }

    bool contains_face(const std::vector<ot::octa_face>& faces, ot::octa_face face)
    {
        return std::ranges::find(faces, face) != faces.end();
    }

    bool contains_vert(const std::vector<ot::octa_vert>& verts, ot::octa_vert vert)
    {
        return std::ranges::find(verts, vert) != verts.end();
    }

    std::pair<ot::octa_vert, ot::octa_vert> canonicalize(ot::octa_vert v1, ot::octa_vert v2)
    {
        return ot::enum_index(v1) < ot::enum_index(v2) ? std::pair(v1, v2) : std::pair(v2, v1);
    }

    ot::cell generate_octa_vert_neighbor_offset(ot::octa_vert vert)
    {
        const ot::cell origin(0, 0, 0);
        std::optional<ot::cell_set> intersection;
        for (ot::octa_edge edge : ot::get_edges_around_octa_vert(vert)) {
            const ot::cell neighbor = ot::get_edge_neighbor(origin, edge);
            ot::cell_set edge_neighbors(neighbor.edge_neighbors());
            intersection = intersection ? intersection->intersect(edge_neighbors) : std::optional<ot::cell_set>(edge_neighbors);
        }

        if (!intersection || intersection->count() != 2) {
            throw std::runtime_error("Could not generate octahedron vertex-neighbor offset");
        }

        return intersection->subtract(ot::cell_set(origin)).cells().front();
    }

    const std::unordered_map<ot::octa_vert, ot::cell>& octa_vert_neighbor_offset_table()
    {
        static const std::unordered_map<ot::octa_vert, ot::cell> table = [] {
            std::unordered_map<ot::octa_vert, ot::cell> output;
            for (ot::octa_vert vert : ot::octa_verts()) {
                output.emplace(vert, generate_octa_vert_neighbor_offset(vert));
            }
            return output;
            }();
        return table;
    }

    std::vector<ot::cell> get_pyramid_neighbors(const ot::cell& index, bool include_pyramid_base = false)
    {
        if (index.type != ot::cell_type::pyramid) {
            throw std::invalid_argument("Expected a pyramid");
        }

        std::vector<ot::cell> output;
        if (include_pyramid_base) {
            output.emplace_back(ot::get_opposite_octa_vert(index.pyramid_tip()), index.column, index.row, index.layer);
        }

        for (ot::octa_face face : ot::octa_faces_around_vert(index.pyramid_tip())) {
            output.push_back(ot::get_octahedron_neighbor(index.column, index.row, index.layer, face));
        }
        return output;
    }

    std::vector<ot::cell> get_tetra_neighbors(
        ot::tetrahedron_orientation orientation,
        int column,
        int row,
        int layer,
        bool include_pyramids = false)
    {
        std::vector<ot::cell> output;
        for (ot::tetra_face_or_vert face : ot::tetra_faces()) {
            const ot::cell octa = ot::get_tetrahedron_neighbor(orientation, column, row, layer, face);
            output.push_back(octa);
            if (include_pyramids) {
                const ot::cell tetra(orientation, column, row, layer);
                for (const ot::cell& pyramid : ot::get_all_pyramids_in_octahedron(octa)) {
                    if (ot::cell_set(pyramid.neighbors()).has_cell(tetra)) {
                        output.push_back(pyramid);
                    }
                }
            }
        }
        return output;
    }

    std::vector<ot::cell> get_octa_neighbors(int column, int row, int layer)
    {
        std::vector<ot::cell> output;
        for (ot::octa_face face : ot::octa_faces()) {
            output.push_back(ot::get_octahedron_neighbor(column, row, layer, face));
        }
        return output;
    }

    bool are_one_cell_apart(const ot::cell& c1, const ot::cell& c2)
    {
        return ot::cell_set(c1.neighbors()).intersect(c2.neighbors()).any();
    }

    std::vector<ot::cell> get_edge_neighbors_impl(const ot::cell& value, bool include_pyramids = false)
    {
        if (value.type == ot::cell_type::tetrahedron) {
            throw std::logic_error("Edge neighbors are not implemented for tetrahedra");
        }

        std::vector<ot::cell> octa_neighbors;
        const std::vector<ot::octa_edge> edges = value.type == ot::cell_type::octahedron
            ? ot::octahedron_edges()
            : ot::pyramid_edges(value.pyramid_tip());

        octa_neighbors.reserve(edges.size());
        for (ot::octa_edge edge : edges) {
            octa_neighbors.push_back(ot::get_edge_neighbor(value, edge));
        }

        if (!include_pyramids) {
            return octa_neighbors;
        }

        std::vector<ot::cell> output = octa_neighbors;
        for (const ot::cell& octa : octa_neighbors) {
            for (const ot::cell& pyramid : ot::get_all_pyramids_in_octahedron(octa)) {
                if (ot::cell_set(pyramid.edge_neighbors()).has_cell(value.to_octahedron()) && are_one_cell_apart(pyramid, value)) {
                    output.push_back(pyramid);
                }
            }
        }
        return output;
    }

    bool nearly_equal(const ot::vector3& lhs, const ot::vector3& rhs)
    {
        return (lhs - rhs).norm() < 0.00000005;
    }

    ot::vector3 origin_octa_vertex_location(ot::octa_vert vert)
    {
        return ot::get_vertex_location(
            ot::get_octahedron_vertex(ot::cell(0, 0, 0), vert));
    }

    bool rotation_maps_directed_edge(
        const ot::matrix3& rotation,
        ot::octa_edge from,
        ot::octa_edge to)
    {
        const auto [from_v0, from_v1] = ot::get_verts_from_octa_edge(from);
        const auto [to_v0, to_v1] = ot::get_verts_from_octa_edge(to);

        const ot::vector3 transformed_0 =
            rotation * origin_octa_vertex_location(from_v0);

        const ot::vector3 transformed_1 =
            rotation * origin_octa_vertex_location(from_v1);

        const ot::vector3 target_0 = origin_octa_vertex_location(to_v0);
        const ot::vector3 target_1 = origin_octa_vertex_location(to_v1);

        return nearly_equal(transformed_0, target_0) &&
            nearly_equal(transformed_1, target_1);
    }

    ot::matrix4 cell_coordinate_matrix_from_space_rotation(const ot::matrix3& rotation)
    {
        const ot::vector3 origin_center = ot::cell(0, 0, 0).center();

        ot::matrix3 cell_to_space;
        cell_to_space.col(0) = ot::cell(1, 0, 0).center() - origin_center;
        cell_to_space.col(1) = ot::cell(0, 1, 0).center() - origin_center;
        cell_to_space.col(2) = ot::cell(0, 0, 1).center() - origin_center;

        const ot::matrix3 cell_rotation =
            cell_to_space.inverse() * rotation * cell_to_space;

        ot::matrix4 output = ot::matrix4::Identity();
        for (int row = 0; row < 3; ++row) {
            for (int column = 0; column < 3; ++column) {
                const double value = cell_rotation(row, column);
                const double rounded = std::round(value);

                if (std::abs(value - rounded) > 0.00000005) {
                    throw std::runtime_error(
                        "Generated octahedron rotation was not integral in cell coordinates");
                }

                output(row, column) = rounded;
            }
        }

        return output;
    }

    constexpr std::size_t octa_edge_count = 12;
    using octahedron_rotation_table =
        std::array<std::array<ot::matrix4, octa_edge_count>, octa_edge_count>;

    ot::matrix4 generate_octahedron_edge_rotation(ot::octa_edge from, ot::octa_edge to)
    {
        if (from == to) {
            return ot::identity_matrix();
        }

        ot::matrix3 octa_axis_to_space;
        octa_axis_to_space.col(0) =
            ot::normalize(origin_octa_vertex_location(ot::octa_vert::top));

        octa_axis_to_space.col(1) =
            ot::normalize(origin_octa_vertex_location(ot::octa_vert::north_east));

        octa_axis_to_space.col(2) =
            ot::normalize(origin_octa_vertex_location(ot::octa_vert::south_east));

        std::array<int, 3> permutation = { 0, 1, 2 };
        do {
            for (int sign_0 : { -1, 1 }) {
                for (int sign_1 : { -1, 1 }) {
                    for (int sign_2 : { -1, 1 }) {
                        const std::array<int, 3> signs = { sign_0, sign_1, sign_2 };

                        ot::matrix3 signed_permutation = ot::matrix3::Zero();
                        for (int i = 0; i < 3; ++i) {
                            signed_permutation(permutation[i], i) =
                                static_cast<double>(signs[i]);
                        }

                        if (signed_permutation.determinant() < 0.0) {
                            continue;
                        }

                        const ot::matrix3 rotation =
                            octa_axis_to_space *
                            signed_permutation *
                            octa_axis_to_space.inverse();

                        if (rotation_maps_directed_edge(rotation, from, to)) {
                            return cell_coordinate_matrix_from_space_rotation(rotation);
                        }
                    }
                }
            }
        } while (std::next_permutation(permutation.begin(), permutation.end()));

        throw std::runtime_error("Could not generate octahedron edge rotation");
    }

    const octahedron_rotation_table& octahedron_edge_rotation_table()
    {
        static const octahedron_rotation_table table = [] {
            octahedron_rotation_table output;
            const std::vector<ot::octa_edge> edges = ot::octahedron_edges();

            if (edges.size() != octa_edge_count) {
                throw std::runtime_error("Unexpected number of octahedron edges");
            }

            for (ot::octa_edge from : edges) {
                const std::size_t from_index =
                    static_cast<std::size_t>(ot::enum_index(from));

                if (from_index >= octa_edge_count) {
                    throw std::runtime_error("Unexpected octahedron edge enum value");
                }

                for (ot::octa_edge to : edges) {
                    const std::size_t to_index =
                        static_cast<std::size_t>(ot::enum_index(to));

                    if (to_index >= octa_edge_count) {
                        throw std::runtime_error("Unexpected octahedron edge enum value");
                    }

                    output[from_index][to_index] =
                        generate_octahedron_edge_rotation(from, to);
                }
            }

            return output;
        }();

        return table;
    }

    int count_octahedron_neighbors_in_set(
        const ot::cell& tetrahedron,
        const ot::cell_set& octahedra)
    {
        if (tetrahedron.type != ot::cell_type::tetrahedron) {
            throw std::invalid_argument("Expected tetrahedron");
        }

        int count = 0;

        for (const ot::cell& neighbor : tetrahedron.neighbors()) {
            if (neighbor.type == ot::cell_type::octahedron &&
                octahedra.has_cell(neighbor)) {
                ++count;
            }
        }

        return count;
    }

} // namespace

const std::array<ot::octa_face, 8>& ot::octa_faces()
{
    static constexpr std::array values = {
        ot::octa_face::north_top,
        ot::octa_face::east_top,
        ot::octa_face::south_top,
        ot::octa_face::west_top,
        ot::octa_face::north_bottom,
        ot::octa_face::east_bottom,
        ot::octa_face::south_bottom,
        ot::octa_face::west_bottom
    };
    return values;
}

const std::array<ot::octa_vert, 6>& ot::octa_verts()
{
    static constexpr std::array values = {
        ot::octa_vert::top,
        ot::octa_vert::north_west,
        ot::octa_vert::north_east,
        ot::octa_vert::south_east,
        ot::octa_vert::south_west,
        ot::octa_vert::bottom
    };
    return values;
}

const std::array<ot::tetra_face_or_vert, 4>& ot::tetra_faces()
{
    static constexpr std::array values = {
        ot::tetra_face_or_vert::north,
        ot::tetra_face_or_vert::south,
        ot::tetra_face_or_vert::east,
        ot::tetra_face_or_vert::west
    };
    return values;
}

const std::array<ot::tetra_face_or_vert, 4>& ot::tetra_verts()
{
    return ot::tetra_faces();
}

std::vector<ot::octa_face> ot::octa_faces_around_vert(ot::octa_vert vert)
{
    return octa_faces_around_vert_table().at(vert);
}

std::vector<ot::octa_vert> ot::get_verts_of_pyramid_base(ot::octa_vert tip)
{
    return octa_base_of_pyramid_table().at(tip);
}

std::vector<ot::octa_edge> ot::octahedron_edges()
{
    std::vector<ot::octa_edge> output;
    output.reserve(octa_edge_to_neighbor_table().size());
    for (const auto& [edge, unused] : octa_edge_to_neighbor_table()) {
        output.push_back(edge);
    }
    return output;
}

std::vector<ot::octa_edge> ot::pyramid_edges(ot::octa_vert tip)
{
    const std::vector<ot::octa_edge> opposite_vert_edges = ot::get_edges_around_octa_vert(ot::get_opposite_octa_vert(tip));
    std::unordered_set<ot::octa_edge> opposite_set(opposite_vert_edges.begin(), opposite_vert_edges.end());

    std::vector<ot::octa_edge> output;
    for (ot::octa_edge edge : ot::octahedron_edges()) {
        if (!opposite_set.contains(edge)) {
            output.push_back(edge);
        }
    }
    return output;
}

ot::cell ot::get_edge_neighbor(const ot::cell& value, ot::octa_edge edge)
{
    if (value.type == ot::cell_type::tetrahedron) {
        throw std::invalid_argument("GetEdgeNeighbor expected an octahedron or pyramid");
    }
    return value + octa_edge_to_neighbor_table().at(edge);
}

std::vector<ot::cell> ot::edge_neighbors(const ot::cell& value, bool include_pyramids)
{
    // The old C# code cached origin-relative answers for the include_pyramids path. The direct implementation is
    // clearer and still cheap for the small local neighborhoods involved here.
    return get_edge_neighbors_impl(value, include_pyramids);
}

ot::cell ot::get_neighbor(const ot::cell& value, ot::octa_face face)
{
    if (value.type == ot::cell_type::tetrahedron) {
        throw std::invalid_argument("GetNeighbor expected an octahedron or pyramid");
    }
    if (value.type == ot::cell_type::pyramid && !contains_face(ot::octa_faces_around_vert(value.pyramid_tip()), face)) {
        throw std::invalid_argument("This pyramid does not contain that face");
    }
    return ot::get_octahedron_neighbor(value.column, value.row, value.layer, face);
}

ot::cell ot::get_neighbor(const ot::cell& value, ot::tetra_face_or_vert face)
{
    if (value.type != ot::cell_type::tetrahedron) {
        throw std::invalid_argument("GetNeighbor expected a tetrahedron");
    }
    return ot::get_tetrahedron_neighbor(value.tetra_orientation(), value.column, value.row, value.layer, face);
}

std::pair<ot::octa_vert, ot::octa_vert> ot::get_verts_from_octa_edge(ot::octa_edge edge)
{
    return octa_edge_to_verts_table().at(edge);
}

ot::octa_edge ot::get_octa_edge_from_verts(ot::octa_vert v1, ot::octa_vert v2)
{
    static const std::map<std::pair<ot::octa_vert, ot::octa_vert>, ot::octa_edge> table = [] {
        std::map<std::pair<ot::octa_vert, ot::octa_vert>, ot::octa_edge> output;
        for (const auto& [edge, verts] : octa_edge_to_verts_table()) {
            output[canonicalize(verts.first, verts.second)] = edge;
        }
        return output;
        }();

    return table.at(canonicalize(v1, v2));
}

std::vector<ot::octa_edge> ot::get_edges_around_octa_vert(ot::octa_vert vert)
{
    static const std::unordered_map<ot::octa_vert, std::vector<ot::octa_edge>> table = [] {
        std::unordered_map<ot::octa_vert, std::vector<ot::octa_edge>> output;
        for (ot::octa_vert value : ot::octa_verts()) {
            output[value] = {};
        }
        for (const auto& [edge, verts] : octa_edge_to_verts_table()) {
            output[verts.first].push_back(edge);
            output[verts.second].push_back(edge);
        }
        return output;
        }();

    return table.at(vert);
}

std::vector<ot::octa_vert> ot::pyramid_vertices(ot::octa_vert tip_of_pyramid)
{
    std::vector<ot::octa_vert> output = { tip_of_pyramid };
    const std::vector<ot::octa_vert> base = ot::get_verts_of_pyramid_base(tip_of_pyramid);
    output.insert(output.end(), base.begin(), base.end());
    return output;
}

ot::cell ot::get_octa_vert_neighbor(const ot::cell& value, ot::octa_vert vert)
{
    if (value.type == ot::cell_type::tetrahedron) {
        throw std::invalid_argument("Expected an octahedron or pyramid");
    }
    if (value.type == ot::cell_type::pyramid && !contains_vert(ot::pyramid_vertices(value.pyramid_tip()), vert)) {
        throw std::invalid_argument("This pyramid does not contain that vertex");
    }
    return value + octa_vert_neighbor_offset_table().at(vert);
}

ot::cell ot::get_octahedron_neighbor(int column, int row, int layer, ot::octa_face face)
{
    ot::cell neighbor(column, row, layer);
    neighbor = neighbor + neighbor_offset_table().at(face);
    return neighbor.to_tetrahedron(neighbor_offset_table().at(face).tetra_orientation());
}

ot::cell ot::get_tetrahedron_neighbor(
    ot::tetrahedron_orientation orientation,
    int column,
    int row,
    int layer,
    ot::tetra_face_or_vert face)
{
    ot::cell neighbor(column, row, layer);
    neighbor = neighbor + tetra_neighbor_offset_table().at({ orientation, face });
    return neighbor.to_octahedron();
}

ot::cell ot::get_octahedron_vertex(const ot::cell& index, ot::octa_vert vert)
{
    switch (vert) {
    case ot::octa_vert::bottom:
        return ot::cell(index.column + 1, index.row + 1, index.layer - 2);
    case ot::octa_vert::north_west:
        return ot::cell(index.column, index.row + 1, index.layer - 1);
    case ot::octa_vert::north_east:
        return ot::cell(index.column + 1, index.row + 1, index.layer - 1);
    case ot::octa_vert::south_east:
        return ot::cell(index.column + 1, index.row, index.layer - 1);
    case ot::octa_vert::south_west:
        return ot::cell(index.column, index.row, index.layer - 1);
    case ot::octa_vert::top:
        return ot::cell(index.column, index.row, index.layer);
    }
    throw std::logic_error("Invalid octahedron vertex");
}

ot::cell ot::get_tetrahedron_vertex(const ot::cell& tetra_index, ot::tetra_face_or_vert vert)
{
    const int column = tetra_index.column;
    const int row = tetra_index.row;
    const int layer = tetra_index.layer;

    if (tetra_index.tetra_orientation() == ot::tetrahedron_orientation::horizontal) {
        switch (vert) {
        case ot::tetra_face_or_vert::north:
            return ot::cell(column, row + 1, layer - 1);
        case ot::tetra_face_or_vert::east:
            return ot::cell(column, row, layer);
        case ot::tetra_face_or_vert::south:
            return ot::cell(column, row, layer - 1);
        case ot::tetra_face_or_vert::west:
            return ot::cell(column - 1, row, layer);
        }
    }
    else {
        switch (vert) {
        case ot::tetra_face_or_vert::north:
            return ot::cell(column, row, layer);
        case ot::tetra_face_or_vert::east:
            return ot::cell(column + 1, row, layer - 1);
        case ot::tetra_face_or_vert::south:
            return ot::cell(column, row - 1, layer);
        case ot::tetra_face_or_vert::west:
            return ot::cell(column, row, layer - 1);
        }
    }
    throw std::logic_error("Illegal tetrahedron vertex");
}

ot::vector3 ot::get_vertex_location(const ot::cell& vert)
{
    return ot::vector3(
        static_cast<double>(vert.column) + 0.5 * static_cast<double>(vert.layer),
        static_cast<double>(vert.row) + 0.5 * static_cast<double>(vert.layer),
        half_root_2 + half_root_2 * static_cast<double>(vert.layer));
}

std::vector<ot::octa_vert> ot::get_verts_for_octa_face(ot::octa_face face, bool use_clockwise_winding_order)
{
    std::vector<ot::octa_vert> verts;
    switch (face) {
    case ot::octa_face::north_top:
        verts = { ot::octa_vert::top, ot::octa_vert::north_west, ot::octa_vert::north_east };
        break;
    case ot::octa_face::east_top:
        verts = { ot::octa_vert::top, ot::octa_vert::north_east, ot::octa_vert::south_east };
        break;
    case ot::octa_face::south_top:
        verts = { ot::octa_vert::top, ot::octa_vert::south_east, ot::octa_vert::south_west };
        break;
    case ot::octa_face::west_top:
        verts = { ot::octa_vert::top, ot::octa_vert::south_west, ot::octa_vert::north_west };
        break;
    case ot::octa_face::north_bottom:
        verts = { ot::octa_vert::bottom, ot::octa_vert::north_east, ot::octa_vert::north_west };
        break;
    case ot::octa_face::east_bottom:
        verts = { ot::octa_vert::bottom, ot::octa_vert::south_east, ot::octa_vert::north_east };
        break;
    case ot::octa_face::south_bottom:
        verts = { ot::octa_vert::bottom, ot::octa_vert::south_west, ot::octa_vert::south_east };
        break;
    case ot::octa_face::west_bottom:
        verts = { ot::octa_vert::bottom, ot::octa_vert::north_west, ot::octa_vert::south_west };
        break;
    }

    if (!use_clockwise_winding_order) {
        std::ranges::reverse(verts);
    }
    return verts;
}

std::vector<ot::tetra_face_or_vert> ot::get_verts_for_tetra_face(
    ot::tetrahedron_orientation orientation,
    ot::tetra_face_or_vert face,
    bool use_clockwise_winding_order)
{
    std::vector<ot::tetra_face_or_vert> verts;
    switch (face) {
    case ot::tetra_face_or_vert::north:
        verts = { ot::tetra_face_or_vert::east, ot::tetra_face_or_vert::west };
        break;
    case ot::tetra_face_or_vert::south:
        verts = { ot::tetra_face_or_vert::west, ot::tetra_face_or_vert::east };
        break;
    case ot::tetra_face_or_vert::east:
        verts = { ot::tetra_face_or_vert::north, ot::tetra_face_or_vert::south };
        break;
    case ot::tetra_face_or_vert::west:
        verts = { ot::tetra_face_or_vert::south, ot::tetra_face_or_vert::north };
        break;
    }

    if (orientation == ot::tetrahedron_orientation::horizontal) {
        verts.insert(verts.begin(), face);
    }
    else {
        std::ranges::reverse(verts);
        verts.insert(verts.begin(), face);
    }

    if (!use_clockwise_winding_order) {
        std::ranges::reverse(verts);
    }
    return verts;
}

ot::octa_vert ot::get_opposite_octa_vert(ot::octa_vert vert)
{
    return opposite_octa_vert_table().at(vert);
}

std::vector<ot::cell> ot::get_all_pyramids_in_octahedron(const ot::cell& octa)
{
    if (octa.type != ot::cell_type::octahedron) {
        throw std::invalid_argument("Expected an octahedron");
    }

    std::vector<ot::cell> output;
    for (ot::octa_vert vert : ot::octa_verts()) {
        output.emplace_back(vert, octa.column, octa.row, octa.layer);
    }
    return output;
}

std::vector<ot::cell> ot::neighbors(const ot::cell& index)
{
    switch (index.type) {
    case ot::cell_type::octahedron:
        return get_octa_neighbors(index.column, index.row, index.layer);
    case ot::cell_type::tetrahedron:
        return get_tetra_neighbors(index.tetra_orientation(), index.column, index.row, index.layer);
    case ot::cell_type::pyramid:
        return get_pyramid_neighbors(index);
    }
    throw std::logic_error("Unknown cell type");
}

std::vector<ot::cell> ot::get_neighbors_with_pyramids(const ot::cell& index, bool include_pyramid_base)
{
    switch (index.type) {
    case ot::cell_type::octahedron:
        return get_octa_neighbors(index.column, index.row, index.layer);
    case ot::cell_type::tetrahedron:
        return get_tetra_neighbors(index.tetra_orientation(), index.column, index.row, index.layer, true);
    case ot::cell_type::pyramid:
        return get_pyramid_neighbors(index, include_pyramid_base);
    }
    throw std::logic_error("Unknown cell type");
}

std::pair<ot::octa_face, ot::octa_face> ot::get_octa_edge(ot::octa_face f1, ot::octa_face f2)
{
    return ot::enum_index(f1) < ot::enum_index(f2) ? std::pair(f1, f2) : std::pair(f2, f1);
}

ot::matrix4 ot::get_cell_rotation_matrix(
    ot::octa_edge column_edge,
    ot::octa_edge row_edge,
    ot::octa_edge layer_edge,
    std::optional<ot::vector3> translation)
{
    const ot::vector3 trans = translation.value_or(ot::vector3::Zero());
    const ot::cell column = octa_edge_to_neighbor_table().at(column_edge);
    const ot::cell row = octa_edge_to_neighbor_table().at(row_edge);
    const ot::cell layer = octa_edge_to_neighbor_table().at(layer_edge);

    ot::matrix4 matrix = ot::matrix4::Identity();
    matrix.col(0) = ot::vector4(column.column, column.row, column.layer, 0.0);
    matrix.col(1) = ot::vector4(row.column, row.row, row.layer, 0.0);
    matrix.col(2) = ot::vector4(layer.column, layer.row, layer.layer, 0.0);
    matrix.col(3) = ot::vector4(trans.x(), trans.y(), trans.z(), 1.0);
    return matrix;
}

ot::matrix4 ot::identity_matrix()
{
    return ot::get_cell_rotation_matrix(ot::octa_edge::east, ot::octa_edge::north, ot::octa_edge::north_east_top);
}

ot::matrix4 ot::get_cell_rotation_matrix(ot::slot edge, int order)
{
    static std::unordered_map<std::pair<ot::slot, int>, ot::matrix4, slot_order_hash> cache;
    const std::pair key(edge, order);
    if (!cache.contains(key)) {
        cache.emplace(key, ot::generate_cell_rotation_matrix(edge, order));
    }
    return cache.at(key);
}

ot::matrix4 ot::generate_cell_rotation_matrix(ot::slot edge, int order)
{
    const bool need_translation = order % 2 == 0;
    switch (edge) {
    case ot::slot::west:
        return ot::get_cell_rotation_matrix(ot::octa_edge::east, ot::octa_edge::north, ot::octa_edge::north_east_top);
    case ot::slot::north:
        return ot::get_cell_rotation_matrix(
            ot::octa_edge::south,
            ot::octa_edge::east,
            ot::octa_edge::south_east_top,
            need_translation ? std::optional(ot::vector3(0.0, -1.0, 0.0)) : std::nullopt);
    case ot::slot::north_west_top:
        return ot::get_cell_rotation_matrix(
            ot::octa_edge::south_east_bottom,
            ot::octa_edge::north_west_bottom,
            ot::octa_edge::north_east_bottom,
            need_translation ? std::optional(ot::vector3(0.0, 0.0, -1.0)) : std::nullopt);
    case ot::slot::north_east_top:
        return ot::get_cell_rotation_matrix(
            ot::octa_edge::south_west_bottom,
            ot::octa_edge::north_east_bottom,
            ot::octa_edge::south_east_bottom,
            need_translation ? std::optional(ot::vector3(0.0, 0.0, -1.0)) : std::nullopt);
    case ot::slot::south_west_top:
        return ot::get_cell_rotation_matrix(ot::octa_edge::north_east_bottom, ot::octa_edge::north_east_top, ot::octa_edge::east);
    case ot::slot::south_east_top:
        return ot::get_cell_rotation_matrix(
            ot::octa_edge::north_west_bottom,
            ot::octa_edge::north_west_top,
            ot::octa_edge::north,
            need_translation ? std::optional(ot::vector3(-1.0, 0.0, 0.0)) : std::nullopt);
    }
    throw std::logic_error("Invalid slot");
}

ot::matrix4 ot::get_flip_matrix(bool lengthwise, int order)
{
    static std::unordered_map<std::pair<bool, int>, ot::matrix4, bool_order_hash> cache;
    const std::pair key(lengthwise, order);
    if (!cache.contains(key)) {
        cache.emplace(key, ot::generate_flip_matrix(lengthwise, order));
    }
    return cache.at(key);
}

ot::matrix4 ot::generate_flip_matrix(bool lengthwise, int order)
{
    const bool need_translation = order % 2 == 0;
    if (lengthwise) {
        return ot::get_cell_rotation_matrix(
            ot::octa_edge::west,
            ot::octa_edge::north,
            ot::octa_edge::north_west_bottom,
            need_translation ? std::optional(ot::vector3(-1.0, 0.0, 0.0)) : std::nullopt);
    }

    return ot::get_cell_rotation_matrix(
        ot::octa_edge::east,
        ot::octa_edge::south,
        ot::octa_edge::south_east_bottom,
        need_translation ? std::optional(ot::vector3(0.0, -1.0, 0.0)) : std::nullopt);
}

std::vector<ot::cell_face> ot::get_tetra_faces(const ot::cell& value)
{
    std::vector<ot::cell_face> output;
    for (ot::tetra_face_or_vert face : ot::tetra_faces()) {
        output.emplace_back(value, face);
    }
    return output;
}

std::vector<ot::cell_face> ot::get_pyramid_faces(const ot::cell& value)
{
    std::vector<ot::cell_face> output;
    for (ot::octa_face face : ot::octa_faces_around_vert(value.pyramid_tip())) {
        output.emplace_back(value, face);
    }
    output.emplace_back(value);
    return output;
}

std::vector<ot::cell_face> ot::get_octa_faces(const ot::cell& value)
{
    std::vector<ot::cell_face> output;
    for (ot::octa_face face : ot::octa_faces()) {
        output.emplace_back(value, face);
    }
    return output;
}

std::vector<ot::cell> ot::explode(const std::vector<ot::cell>& input)
{
    std::unordered_set<ot::cell, ot::cell_hash> set;
    for (const ot::cell& value : input) {
        for (const ot::cell& neighbor : value.neighbors()) {
            set.insert(neighbor);
        }
    }
    return { set.begin(), set.end() };
}

ot::matrix4 ot::translate_octahedra(ot::octa_edge direction, int distance)
{
    const ot::cell offset =
        ot::get_edge_neighbor(ot::cell(0, 0, 0), direction);

    return ot::get_cell_rotation_matrix(
        ot::octa_edge::east,
        ot::octa_edge::north,
        ot::octa_edge::north_east_top,
        ot::vector3(
            static_cast<double>(distance * offset.column),
            static_cast<double>(distance * offset.row),
            static_cast<double>(distance * offset.layer)));
}

ot::cell_set ot::fill_tetrahedra(const ot::cell_set& cells)
{
    ot::cell_set output = cells;

    const ot::cell_set octahedra(cells.octahedra());

    ot::cell_set candidate_tetrahedra;
    for (const ot::cell& octahedron : octahedra.cells()) {
        for (ot::octa_face face : ot::octa_faces()) {
            candidate_tetrahedra.add_cell(
                ot::get_neighbor(octahedron, face));
        }
    }

    for (const ot::cell& tetrahedron : candidate_tetrahedra.cells()) {
        if (count_octahedron_neighbors_in_set(tetrahedron, octahedra) >= 2) {
            output.add_cell(tetrahedron);
        }
    }

    return output;
}

ot::matrix4 ot::rotate_octahedra(ot::octa_edge from, ot::octa_edge to)
{
    return octahedron_edge_rotation_table()
        [static_cast<std::size_t>(ot::enum_index(from))]
        [static_cast<std::size_t>(ot::enum_index(to))];
}
