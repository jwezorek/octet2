#include "cell.hpp"
#include "cell_set.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

ot::cell::cell() = default;

ot::cell::cell(int column_, int row_, int layer_)
    : column(column_)
    , row(row_)
    , layer(layer_)
    , type(ot::cell_type::octahedron)
{}

ot::cell::cell(ot::tetrahedron_orientation orientation, int column_, int row_, int layer_)
    : column(column_)
    , row(row_)
    , layer(layer_)
    , type(ot::cell_type::tetrahedron)
    , tetra_orientation_(orientation)
{}

ot::cell::cell(ot::octa_vert pyramid_tip, int column_, int row_, int layer_)
    : column(column_)
    , row(row_)
    , layer(layer_)
    , type(ot::cell_type::pyramid)
    , pyramid_tip_(pyramid_tip)
{}

ot::octa_vert ot::cell::pyramid_tip() const
{
    if (type != ot::cell_type::pyramid || !pyramid_tip_) {
        throw std::logic_error("Expected a pyramid");
    }
    return *pyramid_tip_;
}

ot::tetrahedron_orientation ot::cell::tetra_orientation() const
{
    if (type != ot::cell_type::tetrahedron || !tetra_orientation_) {
        throw std::logic_error("Expected a tetrahedron");
    }
    return *tetra_orientation_;
}

ot::cell ot::cell::to_octahedron() const
{
    return ot::cell(column, row, layer);
}

ot::cell ot::cell::to_tetrahedron(ot::tetrahedron_orientation orientation) const
{
    return ot::cell(orientation, column, row, layer);
}

ot::cell ot::cell::to_pyramid(ot::octa_vert vert) const
{
    return ot::cell(vert, column, row, layer);
}

std::vector<ot::cell> ot::cell::neighbors() const
{
    return ot::neighbors(*this);
}

std::vector<ot::cell> ot::cell::neighbors_with_pyramids() const
{
    return ot::get_neighbors_with_pyramids(*this);
}

std::vector<ot::cell> ot::cell::edge_neighbors() const
{
    return ot::edge_neighbors(*this);
}

std::vector<ot::cell> ot::cell::edge_neighbors_with_pyramids() const
{
    return ot::edge_neighbors(*this, true);
}

ot::cell ot::cell::translate(int column_delta, int row_delta, int layer_delta) const
{
    switch (type) {
    case ot::cell_type::octahedron:
        return ot::cell(column + column_delta, row + row_delta, layer + layer_delta);
    case ot::cell_type::tetrahedron:
        return ot::cell(tetra_orientation(), column + column_delta, row + row_delta, layer + layer_delta);
    case ot::cell_type::pyramid:
        return ot::cell(pyramid_tip(), column + column_delta, row + row_delta, layer + layer_delta);
    }
    throw std::logic_error("Unknown cell type");
}

ot::cell ot::cell::apply_matrix_to_tetrahedron(const ot::matrix4& matrix) const
{
    std::vector<ot::cell> octahedron_neighbors;
    for (const ot::cell& neighbor : neighbors()) {
        octahedron_neighbors.push_back(neighbor.apply(matrix));
    }

    std::optional<ot::cell_set> intersection;
    for (const ot::cell& octa : octahedron_neighbors) {
        ot::cell_set octa_neighbors(octa.neighbors());
        intersection = intersection ? intersection->intersect(octa_neighbors) : std::optional<ot::cell_set>(octa_neighbors);
    }

    if (intersection && intersection->count() == 1) {
        return intersection->cells().front();
    }

    throw std::runtime_error("Could not transform tetrahedral cell by matrix");
}

ot::cell ot::cell::apply_matrix_to_pyramid(const ot::matrix4& matrix) const
{
    const ot::cell tip_neighbor = ot::get_octa_vert_neighbor(*this, pyramid_tip());
    const ot::cell this_as_octa = to_octahedron();
    const ot::cell rotated_tip_neighbor = tip_neighbor.apply(matrix);
    const ot::cell rotated_this_as_octa = this_as_octa.apply(matrix);

    for (ot::octa_vert vert : ot::octa_verts()) {
        if (ot::get_octa_vert_neighbor(rotated_this_as_octa, vert) == rotated_tip_neighbor) {
            return rotated_this_as_octa.to_pyramid(vert);
        }
    }

    throw std::runtime_error("Could not transform pyramid cell by matrix");
}

ot::cell ot::cell::apply(const ot::matrix4& matrix) const
{
    switch (type) {
    case ot::cell_type::octahedron: {
        const ot::vector4 transformed = matrix * ot::vector4(column, row, layer, 1.0);
        return ot::cell(
            static_cast<int>(std::llround(transformed.x())),
            static_cast<int>(std::llround(transformed.y())),
            static_cast<int>(std::llround(transformed.z())));
    }
    case ot::cell_type::tetrahedron:
        return apply_matrix_to_tetrahedron(matrix);
    case ot::cell_type::pyramid:
        return apply_matrix_to_pyramid(matrix);
    }
    throw std::logic_error("Unknown cell type");
}

std::vector<ot::cell_face> ot::cell::faces() const
{
    switch (type) {
    case ot::cell_type::tetrahedron:
        return ot::get_tetra_faces(*this);
    case ot::cell_type::octahedron:
        return ot::get_octa_faces(*this);
    case ot::cell_type::pyramid:
        return ot::get_pyramid_faces(*this);
    }
    throw std::logic_error("Unknown cell type");
}

std::vector<ot::cell> ot::cell::vertex_indices() const
{
    switch (type) {
    case ot::cell_type::octahedron: {
        std::vector<ot::cell> output;
        for (ot::octa_vert vert : ot::octa_verts()) {
            output.push_back(ot::get_octahedron_vertex(*this, vert));
        }
        return output;
    }
    case ot::cell_type::tetrahedron: {
        std::vector<ot::cell> output;
        for (ot::tetra_face_or_vert vert : ot::tetra_verts()) {
            output.push_back(ot::get_tetrahedron_vertex(*this, vert));
        }
        return output;
    }
    case ot::cell_type::pyramid: {
        std::vector<ot::cell> output;
        for (ot::octa_vert vert : ot::pyramid_vertices(pyramid_tip())) {
            output.push_back(ot::get_octahedron_vertex(*this, vert));
        }
        return output;
    }
    }
    throw std::logic_error("Unknown cell type");
}

ot::vector3 ot::cell::center() const
{
    const ot::cell value = type != ot::cell_type::pyramid ? *this : to_octahedron();
    const std::vector<ot::cell> verts = value.vertex_indices();
    ot::vector3 sum = ot::vector3::Zero();
    for (const ot::cell& vert : verts) {
        sum += ot::get_vertex_location(vert);
    }
    return sum / static_cast<double>(verts.size());
}

ot::cell ot::operator+(const ot::cell& lhs, const ot::cell& rhs)
{
    return ot::cell(lhs.column + rhs.column, lhs.row + rhs.row, lhs.layer + rhs.layer);
}

bool ot::operator==(const ot::cell& lhs, const ot::cell& rhs)
{
    return lhs.column == rhs.column &&
        lhs.row == rhs.row &&
        lhs.layer == rhs.layer &&
        lhs.type == rhs.type &&
        lhs.tetra_orientation_ == rhs.tetra_orientation_ &&
        lhs.pyramid_tip_ == rhs.pyramid_tip_;
}

std::size_t ot::cell_hash::operator()(const ot::cell& value) const noexcept
{
    std::size_t seed = 17;
    combine(seed, std::hash<int>{}(value.column));
    combine(seed, std::hash<int>{}(value.row));
    combine(seed, std::hash<int>{}(value.layer));
    combine(seed, std::hash<int>{}(ot::enum_index(value.type)));
    if (value.type == ot::cell_type::tetrahedron) {
        combine(seed, std::hash<int>{}(ot::enum_index(value.tetra_orientation())));
    }
    if (value.type == ot::cell_type::pyramid) {
        combine(seed, std::hash<int>{}(ot::enum_index(value.pyramid_tip())));
    }
    return seed;
}

void ot::cell_hash::combine(std::size_t& seed, std::size_t value) noexcept
{
    seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
}

ot::cell_face::cell_face(ot::cell value, ot::octa_face face)
    : cell_(value)
    , octa_face_(face)
{
    if (value.type == ot::cell_type::tetrahedron) {
        throw std::invalid_argument("Expected pyramid or octahedron");
    }
}

ot::cell_face::cell_face(ot::cell value, ot::tetra_face_or_vert face)
    : cell_(value)
    , tetra_face_(face)
{
    if (value.type != ot::cell_type::tetrahedron) {
        throw std::invalid_argument("Expected tetrahedron");
    }
}

ot::cell_face::cell_face(ot::cell value)
    : cell_(value)
    , is_square_face_(true)
{
    if (value.type != ot::cell_type::pyramid) {
        throw std::invalid_argument("Expected pyramid");
    }
}

const ot::cell& ot::cell_face::owner_cell() const
{
    return cell_;
}

bool ot::cell_face::is_square_face() const
{
    return is_square_face_;
}

bool ot::cell_face::is_from_tetrahedron() const
{
    return cell_.type == ot::cell_type::tetrahedron;
}

std::vector<ot::cell> ot::cell_face::vertex_indices() const
{
    switch (cell_.type) {
    case ot::cell_type::octahedron: {
        std::vector<ot::cell> output;
        for (ot::octa_vert vert : ot::get_verts_for_octa_face(*octa_face_)) {
            output.push_back(ot::get_octahedron_vertex(cell_, vert));
        }
        return output;
    }
    case ot::cell_type::tetrahedron: {
        std::vector<ot::cell> output;
        for (ot::tetra_face_or_vert vert : ot::get_verts_for_tetra_face(cell_.tetra_orientation(), *tetra_face_)) {
            output.push_back(ot::get_tetrahedron_vertex(cell_, vert));
        }
        return output;
    }
    case ot::cell_type::pyramid: {
        std::vector<ot::cell> output;
        if (is_square_face_) {
            for (ot::octa_vert vert : ot::get_verts_of_pyramid_base(cell_.pyramid_tip())) {
                output.push_back(ot::get_octahedron_vertex(cell_, vert));
            }
        }
        else {
            for (ot::octa_vert vert : ot::get_verts_for_octa_face(*octa_face_)) {
                output.push_back(ot::get_octahedron_vertex(cell_, vert));
            }
        }
        return output;
    }
    }
    throw std::logic_error("Unknown cell type");
}

std::vector<ot::cell_face::edge> ot::cell_face::edges() const
{
    std::vector<ot::cell_face::edge> output;
    const std::size_t n = vertex_indices().size();
    output.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        output.emplace_back(*this, static_cast<int>(i), static_cast<int>((i + 1) % n));
    }
    return output;
}

std::vector<ot::vector3> ot::cell_face::vertices() const
{
    std::vector<ot::vector3> output;
    for (const ot::cell& index : vertex_indices()) {
        output.push_back(ot::get_vertex_location(index));
    }
    return output;
}

bool ot::operator==(const ot::cell_face& lhs, const ot::cell_face& rhs)
{
    return lhs.cell_ == rhs.cell_ &&
        lhs.octa_face_ == rhs.octa_face_ &&
        lhs.tetra_face_ == rhs.tetra_face_ &&
        lhs.is_square_face_ == rhs.is_square_face_;
}

ot::cell_face::edge::edge(ot::cell_face face_, int from_, int to_)
    : face(std::move(face_))
    , from(from_)
    , to(to_)
{}

std::vector<ot::vector3> ot::cell_face::edge::to_vertices() const
{
    const std::vector<ot::vector3> verts = face.vertices();
    return { verts[static_cast<std::size_t>(from)], verts[static_cast<std::size_t>(to)] };
}

ot::vector3 ot::cell_face::edge::from_vertex() const
{
    return to_vertices()[0];
}

ot::vector3 ot::cell_face::edge::to_vertex() const
{
    return to_vertices()[1];
}

std::size_t ot::cell_face_hash::operator()(const ot::cell_face& face) const noexcept
{
    std::size_t seed = 17;
    ot::cell_hash cell_hasher;
    combine(seed, cell_hasher(face.owner_cell()));
    combine(seed, std::hash<bool>{}(face.is_square_face()));
    const auto verts = face.vertex_indices();
    for (const ot::cell& value : verts) {
        combine(seed, cell_hasher(value));
    }
    return seed;
}

void ot::cell_face_hash::combine(std::size_t& seed, std::size_t value) noexcept
{
    seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
}
