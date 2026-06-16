#pragma once

#include "cell.hpp"
#include "cell_set_polyhedron.hpp"
#include "color.hpp"
#include "colored_cell_set.hpp"

#include <string>
#include <unordered_map>
#include <vector>

namespace ot {

    class colored_cell_set_polyhedron {
    public:
        explicit colored_cell_set_polyhedron(float buffer_percent = 1.0f);
        explicit colored_cell_set_polyhedron(const colored_cell_set<ot::color>& cells, float buffer_percent = 1.0f);

        float buffer_percent() const;
        void set_buffer_percent(float value);

        std::vector<ot::color> colors() const;

        void add_cell(const cell& value, ot::color color);
        void write_to_ply_file(const std::string& path) const;

    private:
        void write_to_ply_file_buffered(const std::string& path) const;

        float buffer_percent_ = 1.0f;
        std::unordered_map<ot::color, ot::cell_set_polyhedron, ot::color_hash> cell_sets_;
    };
}
