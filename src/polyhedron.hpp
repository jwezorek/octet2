#pragma once

#include "geometry.hpp"
#include "octet.hpp"
#include "vector_key.hpp"

#include <optional>
#include <string>
#include <vector>

namespace ot {

    class polyhedron {
    public:
        polyhedron();

        void add_face(const std::vector<ot::vector3>& face);

        std::vector<std::vector<int>> faces() const;
        std::vector<ot::vector3> vertices() const;

        void write_to_ply_file(
            const std::string& filename,
            const std::optional<std::string>& comment = std::nullopt) const;

    private:
        int add_vertex(const ot::vector3& vertex);

        vector_dictionary<int> vertex_to_index_;
        std::vector<ot::vector3> vertices_;
        std::vector<std::vector<int>> faces_;
    };

}
