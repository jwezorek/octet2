#include "cell_set.hpp"
#include <stdexcept>

ot::cell_set::cell_set() = default;

ot::cell_set::cell_set(std::optional<ot::cell> value)
{
    if (value) {
        add_cell(*value);
    }
}

ot::cell_set::cell_set(const ot::cell& value)
{
    add_cell(value);
}

ot::cell_set::cell_set(const std::vector<ot::cell>& values)
{
    add_cells(values);
}

std::vector<ot::cell> ot::cell_set::cells() const
{
    return { cells_.begin(), cells_.end() };
}

std::vector<ot::cell> ot::cell_set::octahedra() const
{
    return select_by_type(ot::cell_type::octahedron);
}

std::vector<ot::cell> ot::cell_set::tetrahedra() const
{
    return select_by_type(ot::cell_type::tetrahedron);
}

std::vector<ot::cell> ot::cell_set::pyramids() const
{
    return select_by_type(ot::cell_type::pyramid);
}

std::vector<ot::cell> ot::cell_set::octahedra_and_pyramids() const
{
    std::vector<ot::cell> output;
    for (const ot::cell& value : cells_) {
        if (value.type != ot::cell_type::tetrahedron) {
            output.push_back(value);
        }
    }
    return output;
}

bool ot::cell_set::has_cell(const ot::cell& value) const
{
    return cells_.contains(value);
}

void ot::cell_set::add_cell(const ot::cell& value)
{
    cells_.insert(value);
}

void ot::cell_set::remove_cell(const ot::cell& value)
{
    cells_.erase(value);
}

void ot::cell_set::add_cells(const std::vector<ot::cell>& values)
{
    for (const ot::cell& value : values) {
        add_cell(value);
    }
}

void ot::cell_set::add_cells(const ot::cell_set& values)
{
    for (const ot::cell& value : values.cells_) {
        add_cell(value);
    }
}

ot::cell_set ot::cell_set::intersect(const ot::cell_set& other) const
{
    ot::cell_set output;
    for (const ot::cell& value : cells_) {
        if (other.cells_.contains(value)) {
            output.add_cell(value);
        }
    }
    return output;
}

ot::cell_set ot::cell_set::intersect(const std::vector<ot::cell>& values) const
{
    return intersect(ot::cell_set(values));
}

ot::cell_set ot::cell_set::unite(const ot::cell_set& other) const
{
    ot::cell_set output(*this);
    output.add_cells(other);
    return output;
}

ot::cell_set ot::cell_set::unite(const std::vector<ot::cell>& values) const
{
    return unite(ot::cell_set(values));
}

ot::cell_set ot::cell_set::subtract(const ot::cell_set& other) const
{
    ot::cell_set output;
    for (const ot::cell& value : cells_) {
        if (!other.cells_.contains(value)) {
            output.add_cell(value);
        }
    }
    return output;
}

ot::cell_set ot::cell_set::translate(int column_delta, int row_delta, int layer_delta) const
{
    ot::cell_set output;
    for (const ot::cell& value : cells_) {
        output.add_cell(value.translate(column_delta, row_delta, layer_delta));
    }
    return output;
}

ot::cell_set ot::cell_set::apply(const ot::matrix4& matrix) const
{
    ot::cell_set output;
    for (const ot::cell& value : cells_) {
        output.add_cell(value.apply(matrix));
    }
    return output;
}

int ot::cell_set::pyramid_weight() const
{
    return 2 * static_cast<int>(octahedra().size()) + static_cast<int>(pyramids().size());
}

int ot::cell_set::weight() const
{
    return pyramid_weight() + static_cast<int>(tetrahedra().size());
}

std::size_t ot::cell_set::count() const
{
    return cells_.size();
}

bool ot::cell_set::any() const
{
    return !cells_.empty();
}

std::optional<ot::cell> ot::cell_set::contents_of_octahedral_cell(const ot::cell& index) const
{
    if (index.type != ot::cell_type::octahedron) {
        throw std::invalid_argument("Expected an octahedron");
    }

    if (cells_.contains(index)) {
        return index;
    }

    for (ot::octa_vert vert : ot::octa_verts()) {
        const ot::cell pyramid = index.to_pyramid(vert);
        if (cells_.contains(pyramid)) {
            return pyramid;
        }
    }

    return std::nullopt;
}

bool ot::cell_set::can_add(const ot::cell_set& values) const
{
    const ot::cell_set tetra_intersection =
        ot::cell_set(tetrahedra()).intersect(ot::cell_set(values.tetrahedra()));

    if (tetra_intersection.any()) {
        return false;
    }

    std::vector<ot::cell> this_octa_cells;
    for (const ot::cell& value : octahedra_and_pyramids()) {
        this_octa_cells.push_back(value.to_octahedron());
    }

    std::vector<ot::cell> other_octa_cells;
    for (const ot::cell& value : values.octahedra_and_pyramids()) {
        other_octa_cells.push_back(value.to_octahedron());
    }

    const ot::cell_set intersection =
        ot::cell_set(this_octa_cells).intersect(ot::cell_set(other_octa_cells));

    if (!intersection.any()) {
        return true;
    }

    for (const ot::cell& octa : intersection.cells()) {
        const std::optional<ot::cell> c1 = contents_of_octahedral_cell(octa);
        const std::optional<ot::cell> c2 = values.contents_of_octahedral_cell(octa);

        if (!c1 || !c2) {
            return false;
        }

        if (c1->type == ot::cell_type::octahedron || c2->type == ot::cell_type::octahedron) {
            return false;
        }

        if (c1->pyramid_tip() != ot::get_opposite_octa_vert(c2->pyramid_tip())) {
            return false;
        }
    }

    return true;
}

std::vector<ot::cell> ot::cell_set::select_by_type(ot::cell_type type) const
{
    std::vector<ot::cell> output;
    for (const ot::cell& value : cells_) {
        if (value.type == type) {
            output.push_back(value);
        }
    }
    return output;
}