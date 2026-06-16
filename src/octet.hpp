#pragma once

#include "geometry.hpp"

#include <array>
#include <compare>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace ot {

    enum class cell_type {
        octahedron,
        tetrahedron,
        pyramid
    };

    enum class tetrahedron_orientation {
        horizontal,
        vertical
    };

    enum class octa_face {
        north_top,
        east_top,
        south_top,
        west_top,
        north_bottom,
        east_bottom,
        south_bottom,
        west_bottom
    };

    enum class octa_vert {
        top,
        north_west,
        north_east,
        south_east,
        south_west,
        bottom
    };

    enum class octa_edge {
        north_west_top,
        north_east_top,
        south_east_top,
        south_west_top,
        north_west_bottom,
        north_east_bottom,
        south_east_bottom,
        south_west_bottom,
        north,
        east,
        south,
        west
    };

    enum class tetra_face_or_vert {
        north,
        east,
        south,
        west
    };

    enum class slot {
        north,
        west,
        north_east_top,
        north_west_top,
        south_west_top,
        south_east_top
    };

    template <typename E>
    constexpr auto enum_index(E value) noexcept
    {
        return static_cast<std::underlying_type_t<E>>(value);
    }

    class cell;
    class cell_face;
    class cell_set;

    const std::array<octa_face, 8>& octa_faces();
    const std::array<octa_vert, 6>& octa_verts();
    const std::array<tetra_face_or_vert, 4>& tetra_faces();
    const std::array<tetra_face_or_vert, 4>& tetra_verts();
    std::vector<octa_face> octa_faces_around_vert(octa_vert vert);
    std::vector<octa_vert> get_verts_of_pyramid_base(octa_vert tip);
    std::vector<octa_edge> octahedron_edges();
    std::vector<octa_edge> pyramid_edges(octa_vert tip);
    cell get_edge_neighbor(const cell& value, octa_edge edge);
    std::vector<cell> edge_neighbors(const cell& value, bool include_pyramids = false);
    cell get_neighbor(const cell& value, octa_face face);
    cell get_neighbor(const cell& value, tetra_face_or_vert face);
    std::pair<octa_vert, octa_vert> get_verts_from_octa_edge(octa_edge edge);
    octa_edge get_octa_edge_from_verts(octa_vert v1, octa_vert v2);
    std::vector<octa_edge> get_edges_around_octa_vert(octa_vert vert);
    std::vector<octa_vert> pyramid_vertices(octa_vert tip_of_pyramid);
    cell get_octa_vert_neighbor(const cell& value, octa_vert vert);
    cell get_octahedron_neighbor(int column, int row, int layer, octa_face face);
    cell get_tetrahedron_neighbor(
        tetrahedron_orientation orientation,
        int column,
        int row,
        int layer,
        tetra_face_or_vert face);
    cell get_octahedron_vertex(const cell& index, octa_vert vert);
    cell get_tetrahedron_vertex(const cell& tetra_index, tetra_face_or_vert vert);
    ot::vector3 get_vertex_location(const cell& vert);
    std::vector<octa_vert> get_verts_for_octa_face(
        octa_face face,
        bool use_clockwise_winding_order = true);
    std::vector<tetra_face_or_vert> get_verts_for_tetra_face(
        tetrahedron_orientation orientation,
        tetra_face_or_vert face,
        bool use_clockwise_winding_order = true);
    octa_vert get_opposite_octa_vert(octa_vert vert);
    std::vector<cell> get_all_pyramids_in_octahedron(const cell& octa);
    std::vector<cell> neighbors(const cell& index);
    std::vector<cell> get_neighbors_with_pyramids(
        const cell& index,
        bool include_pyramid_base = false);
    std::pair<octa_face, octa_face> get_octa_edge(octa_face f1, octa_face f2);
    ot::matrix4 get_cell_rotation_matrix(
        octa_edge column_edge,
        octa_edge row_edge,
        octa_edge layer_edge,
        std::optional<ot::vector3> translation = std::nullopt);
    ot::matrix4 identity_matrix();
    ot::matrix4 get_cell_rotation_matrix(slot edge, int order);
    ot::matrix4 generate_cell_rotation_matrix(slot edge, int order);
    ot::matrix4 get_flip_matrix(bool lengthwise, int order);
    ot::matrix4 generate_flip_matrix(bool lengthwise, int order);
    std::vector<cell_face> get_tetra_faces(const cell& value);
    std::vector<cell_face> get_pyramid_faces(const cell& value);
    std::vector<cell_face> get_octa_faces(const cell& value);
    std::vector<cell> explode(const std::vector<cell>& input);
    matrix4 rotate_octahedra(ot::octa_edge from, ot::octa_edge to);
    matrix4 translate_octahedra(ot::octa_edge direction, int distance);
}



