#include "colored_cell_set_polyhedron.hpp"

#include "colored_polyhedron.hpp"

ot::colored_cell_set_polyhedron::colored_cell_set_polyhedron(float buffer_percent)
    : buffer_percent_(buffer_percent)
{}

ot::colored_cell_set_polyhedron::colored_cell_set_polyhedron(
    const ot::colored_cell_set<ot::color>& cells,
    float buffer_percent)
    : colored_cell_set_polyhedron(buffer_percent)
{
    for (const auto& [value, color] : cells.colored_cells()) {
        add_cell(value, color);
    }
}

float ot::colored_cell_set_polyhedron::buffer_percent() const
{
    return buffer_percent_;
}

void ot::colored_cell_set_polyhedron::set_buffer_percent(float value)
{
    buffer_percent_ = value;
    for (auto& [color, cell_set] : cell_sets_) {
        cell_set.set_buffer_percent(value);
    }
}

std::vector<ot::color> ot::colored_cell_set_polyhedron::colors() const
{
    std::vector<ot::color> output;
    output.reserve(cell_sets_.size());
    for (const auto& [color, unused] : cell_sets_) {
        output.push_back(color);
    }
    return output;
}

void ot::colored_cell_set_polyhedron::add_cell(const ot::cell& value, ot::color color)
{
    auto [iter, inserted] = cell_sets_.try_emplace(color, buffer_percent_);
    iter->second.add_cell(value);
}

void ot::colored_cell_set_polyhedron::write_to_ply_file(const std::string& path) const
{
    if (buffer_percent_ < 1.0f) {
        write_to_ply_file_buffered(path);
        return;
    }

    ot::colored_polyhedron poly;
    for (const auto& [color, cell_set] : cell_sets_) {
        for (const ot::cell_face& face : cell_set.faces()) {
            poly.add_face(face.vertices(), color);
        }
    }
    poly.write_to_ply_file(path);
}

void ot::colored_cell_set_polyhedron::write_to_ply_file_buffered(const std::string& path) const
{
    ot::colored_polyhedron poly;
    for (const auto& [color, cell_set] : cell_sets_) {
        for (const std::vector<ot::vector3>& face : cell_set.buffered_faces()) {
            poly.add_face(face, color);
        }
    }
    poly.write_to_ply_file(path);
}
