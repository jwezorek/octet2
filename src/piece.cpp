#include "piece.hpp"
#include "octet.hpp"

namespace {



}

/*------------------------------------------------------------------------------------------------*/

ot::cell_set ot::create_piece(int len, int order) {
    cell_set piece;

    // "left" edge
    cell prev{ 0,0,0 };
    for (int i = 0; i < len; ++i) {
        auto next_cell = get_edge_neighbor(prev, octa_edge::east);
        piece.add_cell(next_cell.to_pyramid(octa_vert::top));
        prev = next_cell;
    }

    // "right" edge
    prev = cell{ 0,0,0 };
    for (int i = 0; i < len; ++i) {
        auto next_cell = get_edge_neighbor(prev, octa_edge::north_west_top);
        piece.add_cell(next_cell.to_pyramid(octa_vert::north_east));
        prev = next_cell;
    }

    // "fill"
    cell prev_row_start = get_edge_neighbor(cell{ 0,0,0 }, octa_edge::east);
    for (int row = 0; row < 3; ++row) {
        auto c = get_edge_neighbor(prev_row_start, octa_edge::north_west_top);
        prev_row_start = c;
        auto row_len = len - row - 1;
        for (int i = 0; i < row_len; ++i) {
            piece.add_cell(c);
            c = get_edge_neighbor(c, octa_edge::east);
        }
    }

    return fill_tetrahedra(piece);
}