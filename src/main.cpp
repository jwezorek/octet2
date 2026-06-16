#include "cell_set.hpp"
#include "octet.hpp"
#include "cell_set_polyhedron.hpp"
#include "colored_cell_set_polyhedron.hpp"
#include "colored_cell_set.hpp"
#include <print>
#include <stdexcept>

namespace {

    ot::cell_set make_equilateral_triangle(int side_length, bool include_tets = true) {

        ot::cell_set output;
        if (side_length <= 0) {
            return {};
        }

        auto lhs = ot::cell(0, 0, 0);
        auto left_offset = ot::get_edge_neighbor(lhs, ot::octa_edge::south_west_bottom);
        auto rhs = ot::cell(0, 0, 0);
        auto right_offset = ot::get_edge_neighbor(rhs, ot::octa_edge::south_east_bottom);

        int last = side_length - 1;

        if (include_tets) {
            output.add_cell(ot::get_neighbor(ot::cell(0, 0, 0), ot::octa_face::north_top));
        }

        for (auto i = 0; i <= last; ++i) {
            for (auto col = lhs.column; col <= rhs.column; ++col) {
                auto c = ot::cell(col, lhs.row, lhs.layer);
                output.add_cell(c);

                if (include_tets) {
                    //if (i != last || col != lhs.column) {
                        output.add_cell(ot::get_neighbor(c, ot::octa_face::west_bottom));
                    //}

                    //if (i != last || col != rhs.column) {
                        output.add_cell(ot::get_neighbor(c, ot::octa_face::east_bottom));
                    //}

                    if (i != side_length - 1) {
                      output.add_cell(ot::get_neighbor(c, ot::octa_face::south_bottom));
                    }
                }
            }
            lhs = lhs + left_offset;
            rhs = rhs + right_offset;
        }

        return output;
    }

    void write_rotation_test_tableau(const std::string& path)
    {
        constexpr ot::octa_edge upper_north_east = ot::octa_edge::north_east_top;
        constexpr ot::octa_edge lower_south_west = ot::octa_edge::south_west_bottom;

        const ot::cell origin(0, 0, 0);

        const ot::cell red_cell =
            origin.apply(ot::translate_octahedra(upper_north_east, 2));

        const ot::cell yellow_cell =
            origin.apply(ot::translate_octahedra(lower_south_west, 2));

        ot::colored_cell_set<ot::color> tableau;
        tableau.add_cell(red_cell, ot::color(0xff0000));
        tableau.add_cell(yellow_cell, ot::color(0xffff00));

        ot::cell_set test;
        test.add_cell(origin);
        test.add_cell(origin.apply(ot::translate_octahedra(upper_north_east, 1)));

        const ot::matrix4 rotation =
            ot::rotate_octahedra(upper_north_east, lower_south_west);

        const ot::cell_set rotated_test = test.apply(rotation);

        for (const ot::cell& value : rotated_test.cells()) {
            tableau.add_cell(value, ot::color(0x0000ff));
        }

        ot::colored_cell_set_polyhedron poly(tableau);
        poly.write_to_ply_file(path);
    }
}

int main()
{
    const ot::cell_set piece_1 = make_equilateral_triangle(5, false);

    ot::colored_cell_set<ot::color> cells;
    cells.add_cells(piece_1.cells(), ot::color(0xffffff));
    cells.add_cell(ot::cell{ 0,0,0 }, ot::color(0xff0000));

    ot::colored_cell_set_polyhedron poly(cells);
    poly.write_to_ply_file("C:\\test\\octet_colored.ply");

    write_rotation_test_tableau("C:\\test\\rotation_test_tableau.ply");

    std::println("wrote colored model");
}
