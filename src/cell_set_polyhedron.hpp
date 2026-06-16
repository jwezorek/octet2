#pragma once

#include "cell.hpp"
#include "cell_set.hpp"
#include "polyhedron.hpp"

#include <string>
#include <unordered_map>
#include <vector>

namespace ot {

    class cell_set_polyhedron {
    public:
        explicit cell_set_polyhedron(float buffer_percent = 1.0f);
        explicit cell_set_polyhedron(const cell_set& cells, float buffer_percent = 1.0f);

        float buffer_percent() const;
        void set_buffer_percent(float value);

        std::vector<cell_face> faces() const;
        std::vector<std::vector<ot::vector3>> buffered_faces() const;
        std::vector<ot::vector3> vertices() const;

        void add_cell(const cell& value);
        void write_to_ply_file(const std::string& path) const;

    private:
        class stitch {
        public:
            ot::vector3 from;
            ot::vector3 to;
            int from_index = -1;
            int to_index = -1;

            stitch(const ot::vector3& from_, const ot::vector3& to_)
                : from(from_)
                , to(to_)
            {}

            std::vector<ot::vector3> vertices() const
            {
                return { from, to };
            }
        };

        using face_vertex_map = std::unordered_map<cell_face, std::vector<ot::vector3>, cell_face_hash>;
        using edge_pair = std::pair<cell_face::edge, cell_face::edge>;

        polyhedron get_buffered_polyhedron() const;
        void write_to_ply_file_buffered(const std::string& path) const;
        vector_dictionary<std::vector<std::vector<ot::vector3>>> generate_vertex_to_buffered_faces_dictionary(
            const vec_list_dictionary<cell_face>& verts_to_face,
            const face_vertex_map& face_to_buffered_verts) const;
        void stitch_vertices(
            polyhedron& poly,
            vector_dictionary<std::vector<stitch>>& stitches_per_vert,
            vector_dictionary<std::vector<std::vector<ot::vector3>>>& vert_to_faces) const;
        static bool is_face(
            const ot::vector3& v1,
            const ot::vector3& v2,
            const ot::vector3& v3,
            const std::vector<std::vector<ot::vector3>>& faces);
        static std::vector<std::vector<int>> get_loops(const std::vector<stitch>& all_stitches, int n);
        static std::vector<int> get_loop(const std::vector<stitch>& stitches, int n);
        static std::vector<stitch> prune_stitches(const std::vector<stitch>& stitches);
        static std::vector<ot::vector3> fill_in_stitch_indices(std::vector<stitch>& stitches);
        static ot::vector3 center_of(const std::vector<ot::vector3>& points);
        static void stitch_edges(
            polyhedron& poly,
            vector_dictionary<std::vector<stitch>>& stitches,
            const face_vertex_map& face_to_buffered_verts,
            const std::vector<edge_pair>& edges);
        static void add_stitch(vector_dictionary<std::vector<stitch>>& stitches, const ot::vector3& vertex, const stitch& value);
        static void stitch_non_coplanar_edge(
            polyhedron& poly,
            vector_dictionary<std::vector<stitch>>& stitches,
            const cell_face::edge& edge1,
            const std::vector<ot::vector3>& face1,
            const cell_face::edge& edge2,
            const std::vector<ot::vector3>& face2);
        static ot::vector3 previous(const std::vector<ot::vector3>& face, int index);
        static ot::vector3 next(const std::vector<ot::vector3>& face, int index);
        static bool are_coplanar(
            const ot::vector3& a,
            const ot::vector3& b,
            const ot::vector3& c,
            const ot::vector3& d,
            const ot::vector3& e);
        static std::vector<edge_pair> get_edge_pairs(const std::vector<cell_face>& faces);
        static face_vertex_map get_buffered_faces(float octa_buffer_percent, const std::vector<cell_face>& faces);

        float buffer_percent_ = 1.0f;
        vec_list_dictionary<cell_face> vertices_to_cell_face_;
    };
}
