#include "colored_polyhedron.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <utility>

ot::colored_polyhedron::colored_polyhedron()
    : vertex_to_index_(0.00005)
{}

void ot::colored_polyhedron::add_face(const std::vector<ot::vector3>& face, ot::color color)
{
    std::vector<int> indexed_face;
    indexed_face.reserve(face.size());
    for (const ot::vector3& vertex : face) {
        indexed_face.push_back(add_vertex(vertex, color));
    }
    std::ranges::reverse(indexed_face);
    faces_.push_back(std::move(indexed_face));
}

std::vector<std::vector<int>> ot::colored_polyhedron::faces() const
{
    return faces_;
}

std::vector<ot::colored_vertex> ot::colored_polyhedron::vertices() const
{
    return vertices_;
}

void ot::colored_polyhedron::write_to_ply_file(
    const std::string& filename,
    const std::optional<std::string>& comment) const
{
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("Unable to open PLY file for writing: " + filename);
    }

    file << "ply\n";
    file << "format ascii 1.0\n";
    if (comment && !comment->empty()) {
        file << "comment " << *comment << '\n';
    }

    file << "element vertex " << vertices_.size() << '\n';
    file << "property float32 x\n";
    file << "property float32 y\n";
    file << "property float32 z\n";
    file << "property uchar red\n";
    file << "property uchar green\n";
    file << "property uchar blue\n";
    file << "element face " << faces_.size() << '\n';
    file << "property list uint8 int32 vertex_index\n";
    file << "end_header\n";

    for (const ot::colored_vertex& vertex : vertices_) {
        file << vertex.position.x() << ' '
             << vertex.position.y() << ' '
             << vertex.position.z() << ' '
             << static_cast<int>(vertex.color.red) << ' '
             << static_cast<int>(vertex.color.green) << ' '
             << static_cast<int>(vertex.color.blue) << '\n';
    }

    for (const std::vector<int>& face : faces_) {
        file << face.size();
        for (int index : face) {
            file << ' ' << index;
        }
        file << '\n';
    }
}

int ot::colored_polyhedron::add_vertex(const ot::vector3& vertex, ot::color color)
{
    if (!vertex_to_index_.contains(vertex, color)) {
        vertex_to_index_(vertex, color) = static_cast<int>(vertices_.size());
        vertices_.push_back(ot::colored_vertex{ vertex, color });
    }
    return vertex_to_index_.at(vertex, color);
}
