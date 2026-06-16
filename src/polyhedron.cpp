#include "polyhedron.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <utility>

ot::polyhedron::polyhedron()
    : vertex_to_index_(0.00005)
{}

void ot::polyhedron::add_face(const std::vector<ot::vector3>& face)
{
    std::vector<int> indexed_face;
    indexed_face.reserve(face.size());
    for (const ot::vector3& vertex : face) {
        indexed_face.push_back(add_vertex(vertex));
    }
    std::ranges::reverse(indexed_face);
    faces_.push_back(std::move(indexed_face));
}

std::vector<std::vector<int>> ot::polyhedron::faces() const
{
    return faces_;
}

std::vector<ot::vector3> ot::polyhedron::vertices() const
{
    return vertices_;
}

void ot::polyhedron::write_to_ply_file(
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
    file << "element face " << faces_.size() << '\n';
    file << "property list uint8 int32 vertex_index\n";
    file << "end_header\n";

    for (const ot::vector3& vertex : vertices_) {
        file << vertex.x() << ' ' << vertex.y() << ' ' << vertex.z() << '\n';
    }

    for (const std::vector<int>& face : faces_) {
        file << face.size();
        for (int index : face) {
            file << ' ' << index;
        }
        file << '\n';
    }
}

int ot::polyhedron::add_vertex(const ot::vector3& vertex)
{
    if (!vertex_to_index_.contains(vertex)) {
        vertex_to_index_[vertex] = static_cast<int>(vertices_.size());
        vertices_.push_back(vertex);
    }
    return vertex_to_index_.at(vertex);
}
