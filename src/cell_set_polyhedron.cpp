#include "cell_set_polyhedron.hpp"

#include <algorithm>
#include <array>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <unordered_set>

namespace {

    constexpr double epsilon = 0.00005;

    struct int_pair_hash {
        std::size_t operator()(const std::pair<int, int>& value) const noexcept
        {
            std::size_t seed = 17;
            seed ^= std::hash<int>{}(value.first) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
            seed ^= std::hash<int>{}(value.second) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
            return seed;
        }
    };

    std::vector<ot::vector3> face_vertices_from_indices(const ot::cell_face& face)
    {
        std::vector<ot::vector3> output;
        for (const ot::cell& index : face.vertex_indices()) {
            output.push_back(ot::get_vertex_location(index));
        }
        return output;
    }

} // namespace

ot::cell_set_polyhedron::cell_set_polyhedron(float buffer_percent)
    : buffer_percent_(buffer_percent)
    , vertices_to_cell_face_(epsilon, false)
{
}

ot::cell_set_polyhedron::cell_set_polyhedron(const ot::cell_set& cells, float buffer_percent)
    : cell_set_polyhedron(buffer_percent)
{
    for (const cell& value : cells.cells()) {
        add_cell(value);
    }
}

float ot::cell_set_polyhedron::buffer_percent() const
{
    return buffer_percent_;
}

void ot::cell_set_polyhedron::set_buffer_percent(float value)
{
    buffer_percent_ = value;
}

std::vector<ot::cell_face> ot::cell_set_polyhedron::faces() const
{
    return vertices_to_cell_face_.values();
}

std::vector<std::vector<ot::vector3>> ot::cell_set_polyhedron::buffered_faces() const
{
    const polyhedron poly = get_buffered_polyhedron();
    const std::vector<ot::vector3> verts = poly.vertices();
    std::vector<std::vector<ot::vector3>> output;
    for (const std::vector<int>& face : poly.faces()) {
        std::vector<ot::vector3> points;
        points.reserve(face.size());
        for (int index : face) {
            points.push_back(verts[static_cast<std::size_t>(index)]);
        }
        output.push_back(std::move(points));
    }
    return output;
}

std::vector<ot::vector3> ot::cell_set_polyhedron::vertices() const
{
    vector_set unique_vertices(epsilon);
    std::vector<ot::vector3> output;
    for (const cell_face& face : vertices_to_cell_face_.values()) {
        for (const ot::vector3& vertex : face.vertices()) {
            if (!unique_vertices.contains(vertex)) {
                unique_vertices.add(vertex);
                output.push_back(vertex);
            }
        }
    }
    return output;
}

ot::cell_set_polyhedron::face_vertex_map ot::cell_set_polyhedron::get_buffered_faces(
    float octa_buffer_percent,
    const std::vector<cell_face>& faces)
{
    const float tetra_buffer_percent = 2.0f * octa_buffer_percent - 1.0f;
    face_vertex_map output;
    for (const cell_face& face : faces) {
        const ot::vector3 center_of_cell = face.owner_cell().center();
        const std::vector<ot::vector3> face_verts = face.vertices();
        const double scale = face.is_from_tetrahedron()
            ? static_cast<double>(tetra_buffer_percent)
            : static_cast<double>(octa_buffer_percent);

        const Eigen::Affine3d transform = Eigen::Translation3d(center_of_cell) *
            Eigen::Scaling(scale) *
            Eigen::Translation3d(-center_of_cell);

        std::vector<ot::vector3> buffered;
        buffered.reserve(face_verts.size());
        for (const ot::vector3& vertex : face_verts) {
            buffered.push_back(ot::transform_point(vertex, transform.matrix()));
        }
        output.emplace(face, std::move(buffered));
    }
    return output;
}

void ot::cell_set_polyhedron::add_cell(const cell& value)
{
    for (const cell_face& face : value.faces()) {
        const std::vector<ot::vector3> face_verts = face.vertices();
        if (vertices_to_cell_face_.contains(face_verts)) {
            vertices_to_cell_face_.remove(face_verts);
        }
        else {
            vertices_to_cell_face_.insert_or_assign(face_verts, face);
        }
    }
}

ot::polyhedron ot::cell_set_polyhedron::get_buffered_polyhedron() const
{
    polyhedron poly;
    const std::vector<cell_face> current_faces = faces();
    const face_vertex_map face_to_buffered_verts = get_buffered_faces(buffer_percent_, current_faces);
    const std::vector<edge_pair> edges = get_edge_pairs(current_faces);

    for (const auto& [face, vertices] : face_to_buffered_verts) {
        poly.add_face(vertices);
    }

    vector_dictionary<std::vector<stitch>> stitches(epsilon);
    stitch_edges(poly, stitches, face_to_buffered_verts, edges);

    auto faces_per_vertex = generate_vertex_to_buffered_faces_dictionary(vertices_to_cell_face_, face_to_buffered_verts);
    stitch_vertices(poly, stitches, faces_per_vertex);
    return poly;
}

void ot::cell_set_polyhedron::write_to_ply_file_buffered(const std::string& path) const
{
    const polyhedron poly = get_buffered_polyhedron();
    poly.write_to_ply_file(path);
}

ot::vector_dictionary<std::vector<std::vector<ot::vector3>>> ot::cell_set_polyhedron::generate_vertex_to_buffered_faces_dictionary(
    const vec_list_dictionary<cell_face>& verts_to_face,
    const face_vertex_map& face_to_buffered_verts) const
{
    vector_dictionary<std::vector<cell_face>> vert_to_faces(epsilon);
    for (const cell_face& face : verts_to_face.values()) {
        for (const ot::vector3& vertex : face.vertices()) {
            if (!vert_to_faces.contains(vertex)) {
                vert_to_faces[vertex] = { face };
            }
            else {
                vert_to_faces[vertex].push_back(face);
            }
        }
    }

    vector_dictionary<std::vector<std::vector<ot::vector3>>> vert_to_buffered_face(epsilon);
    for (const ot::vector3& vertex : vertices()) {
        std::vector<std::vector<ot::vector3>> buffered_verts;
        for (const cell_face& face : vert_to_faces.at(vertex)) {
            buffered_verts.push_back(face_to_buffered_verts.at(face));
        }
        vert_to_buffered_face[vertex] = std::move(buffered_verts);
    }
    return vert_to_buffered_face;
}

void ot::cell_set_polyhedron::stitch_vertices(
    polyhedron& poly,
    vector_dictionary<std::vector<stitch>>& stitches_per_vert,
    vector_dictionary<std::vector<std::vector<ot::vector3>>>& vert_to_faces) const
{
    for (const ot::vector3& vertex : vertices()) {
        if (!stitches_per_vert.contains(vertex)) {
            continue;
        }

        std::vector<stitch> stitches = stitches_per_vert[vertex];
        std::vector<ot::vector3> stitch_vertices = fill_in_stitch_indices(stitches);
        stitches = prune_stitches(stitches);
        if (stitches.empty()) {
            continue;
        }

        const std::vector<std::vector<int>> loops = get_loops(stitches, static_cast<int>(stitch_vertices.size()));
        for (const std::vector<int>& initial_loop : loops) {
            std::vector<int> loop = initial_loop;
            if (loop.size() == 3) {
                poly.add_face({ stitch_vertices[loop[0]], stitch_vertices[loop[1]], stitch_vertices[loop[2]] });
                continue;
            }

            while (!loop.empty()) {
                bool added_face = false;
                for (std::size_t i = 0; !added_face && i < loop.size(); ++i) {
                    const int previous_index = i > 0 ? loop[i - 1] : loop.back();
                    const int node = loop[i];
                    const int next_index = i + 1 < loop.size() ? loop[i + 1] : loop.front();

                    if (is_face(
                        stitch_vertices[static_cast<std::size_t>(previous_index)],
                        stitch_vertices[static_cast<std::size_t>(node)],
                        stitch_vertices[static_cast<std::size_t>(next_index)],
                        vert_to_faces[vertex])) {
                        poly.add_face({
                            stitch_vertices[static_cast<std::size_t>(previous_index)],
                            stitch_vertices[static_cast<std::size_t>(node)],
                            stitch_vertices[static_cast<std::size_t>(next_index)]
                        });
                        loop.erase(loop.begin() + static_cast<std::ptrdiff_t>(i));
                        added_face = true;
                    }
                }

                if (loop.size() == 3) {
                    poly.add_face({ stitch_vertices[loop[0]], stitch_vertices[loop[1]], stitch_vertices[loop[2]] });
                    loop.clear();
                }

                if (!added_face && !loop.empty()) {
                    throw std::runtime_error("Could not triangulate buffered vertex stitch loop");
                }
            }
        }
    }
}

bool ot::cell_set_polyhedron::is_face(
    const ot::vector3& v1,
    const ot::vector3& v2,
    const ot::vector3& v3,
    const std::vector<std::vector<ot::vector3>>& faces)
{
    for (const std::vector<ot::vector3>& face : faces) {
        if (ot::are_coplanar(v1, v2, v3, center_of(face), epsilon)) {
            return true;
        }
    }
    return false;
}

std::vector<std::vector<int>> ot::cell_set_polyhedron::get_loops(const std::vector<stitch>& all_stitches, int n)
{
    std::vector<std::vector<int>> loops;
    std::vector<stitch> stitches = all_stitches;
    while (!stitches.empty()) {
        const std::vector<int> loop = get_loop(stitches, n);
        loops.push_back(loop);

        std::vector<std::pair<int, int>> pairs;
        for (std::size_t i = 0; i + 1 < loop.size(); ++i) {
            pairs.emplace_back(loop[i], loop[i + 1]);
        }
        pairs.emplace_back(loop.back(), loop.front());
        std::unordered_set<std::pair<int, int>, int_pair_hash> used(pairs.begin(), pairs.end());

        std::vector<stitch> remaining;
        for (const stitch& value : stitches) {
            if (!used.contains({ value.from_index, value.to_index })) {
                remaining.push_back(value);
            }
        }
        stitches = std::move(remaining);
    }
    return loops;
}

std::vector<int> ot::cell_set_polyhedron::get_loop(const std::vector<stitch>& stitches, int n)
{
    std::vector<std::optional<stitch>> stitch_by_from(static_cast<std::size_t>(n));
    for (const stitch& value : stitches) {
        stitch_by_from[static_cast<std::size_t>(value.from_index)] = value;
    }

    const int start = stitches.front().from_index;
    std::vector<int> loop = { start };
    std::optional<stitch> current = stitch_by_from[static_cast<std::size_t>(start)];
    std::unordered_set<int> visited;
    while (current && current->to_index != start) {
        loop.push_back(current->to_index);
        current = stitch_by_from[static_cast<std::size_t>(current->to_index)];
        if (!current) {
            break;
        }
        visited.insert(current->from_index);
        if (visited.contains(current->to_index)) {
            break;
        }
    }

    if (current && current->to_index != start) {
        auto iter = std::ranges::find(loop, current->to_index);
        if (iter != loop.end()) {
            loop.erase(loop.begin(), iter);
        }
    }
    return loop;
}

std::vector<ot::cell_set_polyhedron::stitch> ot::cell_set_polyhedron::prune_stitches(const std::vector<stitch>& stitches)
{
    std::unordered_set<std::pair<int, int>, int_pair_hash> set;
    for (const stitch& value : stitches) {
        set.emplace(value.from_index, value.to_index);
    }

    std::vector<stitch> output;
    for (const stitch& value : stitches) {
        if (!set.contains({ value.to_index, value.from_index })) {
            output.push_back(value);
        }
    }
    return output;
}

std::vector<ot::vector3> ot::cell_set_polyhedron::fill_in_stitch_indices(std::vector<stitch>& stitches)
{
    vector_dictionary<int> vert_to_index(epsilon);
    std::vector<ot::vector3> vertices;
    for (stitch& value : stitches) {
        for (const ot::vector3& vertex : value.vertices()) {
            if (!vert_to_index.contains(vertex)) {
                vert_to_index[vertex] = static_cast<int>(vertices.size());
                vertices.push_back(vertex);
            }
        }
    }

    for (stitch& value : stitches) {
        value.from_index = vert_to_index.at(value.from);
        value.to_index = vert_to_index.at(value.to);
    }
    return vertices;
}

ot::vector3 ot::cell_set_polyhedron::center_of(const std::vector<ot::vector3>& points)
{
    if (points.empty()) {
        throw std::invalid_argument("center_of requires at least one point");
    }

    ot::vector3 sum = ot::vector3::Zero();
    for (const ot::vector3& point : points) {
        sum += point;
    }
    return sum / static_cast<double>(points.size());
}

void ot::cell_set_polyhedron::stitch_edges(
    polyhedron& poly,
    vector_dictionary<std::vector<stitch>>& stitches,
    const face_vertex_map& face_to_buffered_verts,
    const std::vector<edge_pair>& edges)
{
    for (const edge_pair& edge : edges) {
        const cell_face::edge& edge1 = edge.first;
        const cell_face::edge& edge2 = edge.second;

        const std::vector<ot::vector3>& face1_buffered = face_to_buffered_verts.at(edge1.face);
        const std::vector<ot::vector3>& face2_buffered = face_to_buffered_verts.at(edge2.face);

        const ot::vector3 a = face1_buffered[static_cast<std::size_t>(edge1.from)];
        const ot::vector3 b = face2_buffered[static_cast<std::size_t>(edge2.to)];

        if ((a - b).norm() < epsilon) {
            continue;
        }

        const ot::vector3 c = face2_buffered[static_cast<std::size_t>(edge2.from)];
        const ot::vector3 d = face1_buffered[static_cast<std::size_t>(edge1.to)];
        const ot::vector3 center1 = center_of(face1_buffered);
        const ot::vector3 center2 = center_of(face2_buffered);

        if (are_coplanar(a, b, c, d, center1) || are_coplanar(a, b, c, d, center2)) {
            poly.add_face({ a, b, d });
            poly.add_face({ b, c, d });
            add_stitch(stitches, edge1.from_vertex(), stitch(b, a));
            add_stitch(stitches, edge1.to_vertex(), stitch(d, c));
        }
        else {
            stitch_non_coplanar_edge(poly, stitches, edge1, face1_buffered, edge2, face2_buffered);
        }
    }
}

void ot::cell_set_polyhedron::add_stitch(
    vector_dictionary<std::vector<stitch>>& stitches,
    const ot::vector3& vertex,
    const stitch& value)
{
    if (!stitches.contains(vertex)) {
        stitches[vertex] = { value };
    }
    else {
        stitches[vertex].push_back(value);
    }
}

void ot::cell_set_polyhedron::stitch_non_coplanar_edge(
    polyhedron& poly,
    vector_dictionary<std::vector<stitch>>& stitches,
    const cell_face::edge& edge1,
    const std::vector<ot::vector3>& face1,
    const cell_face::edge& edge2,
    const std::vector<ot::vector3>& face2)
{
    const ot::vector3 a = face1[static_cast<std::size_t>(edge1.from)];
    const ot::vector3 b = face2[static_cast<std::size_t>(edge2.to)];
    const ot::vector3 c = face2[static_cast<std::size_t>(edge2.from)];
    const ot::vector3 d = face1[static_cast<std::size_t>(edge1.to)];
    const ot::vector3 before_a = previous(face1, edge1.from);
    const ot::vector3 after_d = next(face1, edge1.to);
    const ot::vector3 after_b = next(face2, edge2.to);
    const ot::vector3 before_c = previous(face2, edge2.from);

    const std::optional<ot::vector3> intersect1 = ot::line_line_intersect(before_a, a, b, after_b, epsilon);
    const std::optional<ot::vector3> intersect2 = ot::line_line_intersect(d, after_d, before_c, c, epsilon);
    if (!intersect1 || !intersect2) {
        throw std::runtime_error("Could not intersect non-coplanar edge stitch lines");
    }

    const ot::vector3 e = *intersect1;
    const ot::vector3 f = *intersect2;

    add_stitch(stitches, edge1.from_vertex(), stitch(b, e));
    add_stitch(stitches, edge1.from_vertex(), stitch(e, a));
    add_stitch(stitches, edge1.to_vertex(), stitch(d, f));
    add_stitch(stitches, edge1.to_vertex(), stitch(f, c));

    for (const std::vector<ot::vector3>& face : std::vector<std::vector<ot::vector3>> {
        { a, e, d },
        { e, f, d },
        { e, b, c },
        { c, f, e }
    }) {
        poly.add_face(face);
    }
}

ot::vector3 ot::cell_set_polyhedron::previous(const std::vector<ot::vector3>& face, int index)
{
    return index > 0 ? face[static_cast<std::size_t>(index - 1)] : face.back();
}

ot::vector3 ot::cell_set_polyhedron::next(const std::vector<ot::vector3>& face, int index)
{
    return index + 1 < static_cast<int>(face.size()) ? face[static_cast<std::size_t>(index + 1)] : face.front();
}

bool ot::cell_set_polyhedron::are_coplanar(
    const ot::vector3& a,
    const ot::vector3& b,
    const ot::vector3& c,
    const ot::vector3& d,
    const ot::vector3& e)
{
    const std::array<ot::vector3, 5> values = { a, b, c, d, e };
    for (std::size_t omitted = 0; omitted < values.size(); ++omitted) {
        std::vector<ot::vector3> combo;
        combo.reserve(4);
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (i != omitted) {
                combo.push_back(values[i]);
            }
        }
        if (!ot::are_coplanar(combo[0], combo[1], combo[2], combo[3], epsilon)) {
            return false;
        }
    }
    return true;
}

std::vector<ot::cell_set_polyhedron::edge_pair> ot::cell_set_polyhedron::get_edge_pairs(const std::vector<cell_face>& faces)
{
    vec_list_dictionary<std::vector<cell_face::edge>> vertices_to_edges(epsilon, false);
    for (const cell_face& face : faces) {
        for (const cell_face::edge& edge : face.edges()) {
            const std::vector<ot::vector3> verts = edge.to_vertices();
            if (vertices_to_edges.contains(verts)) {
                vertices_to_edges[verts].push_back(edge);
            }
            else {
                vertices_to_edges[verts] = { edge };
            }
        }
    }

    std::vector<edge_pair> output;
    for (const std::vector<cell_face::edge>& list : vertices_to_edges.values()) {
        if (list.size() != 2) {
            throw std::runtime_error("Non-manifold edge while building cell-set polyhedron");
        }
        output.emplace_back(list[0], list[1]);
    }
    return output;
}

void ot::cell_set_polyhedron::write_to_ply_file(const std::string& path) const
{
    if (buffer_percent_ < 1.0f) {
        write_to_ply_file_buffered(path);
        return;
    }

    polyhedron poly;
    for (const cell_face& face : faces()) {
        poly.add_face(face_vertices_from_indices(face));
    }
    poly.write_to_ply_file(path);
}

