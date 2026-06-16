#pragma once

#include "color.hpp"
#include "colored_vector_dictionary.hpp"
#include "geometry.hpp"

#include <optional>
#include <string>
#include <vector>

namespace ot {

    struct colored_vertex {
        ot::vector3 position;
        ot::color color;
    };

    class colored_polyhedron {
    public:
        colored_polyhedron();

        void add_face(const std::vector<ot::vector3>& face, ot::color color);

        std::vector<std::vector<int>> faces() const;
        std::vector<ot::colored_vertex> vertices() const;

        void write_to_ply_file(
            const std::string& filename,
            const std::optional<std::string>& comment = std::nullopt) const;

    private:
        int add_vertex(const ot::vector3& vertex, ot::color color);

        colored_vector_dictionary<int> vertex_to_index_;
        std::vector<ot::colored_vertex> vertices_;
        std::vector<std::vector<int>> faces_;
    };
}
