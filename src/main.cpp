#include "cell_set.hpp"
#include "octet.hpp"
#include "cell_set_polyhedron.hpp"
#include "colored_cell_set_polyhedron.hpp"
#include "colored_cell_set.hpp"
#include "piece.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <print>
#include <stdexcept>
#include <string>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

namespace {

    constexpr std::size_t triangle_count = 4;

    // Used only for the diagnostic one-triangle red/white outputs.
    // This shrinks cells slightly so the red cells are easier to see.
    constexpr float diagnostic_buffer_percent = 1.0f;

    ot::cell global_center_octahedron()
    {
        return ot::cell(0, 0, 0);
    }

    struct cuboct_triangle_model {
        std::array<ot::cell_set, triangle_count> full_triangles;
        std::array<ot::cell_set, triangle_count> octahedron_triangles;

        ot::colored_cell_set<ot::color> full_compound;
        ot::colored_cell_set<ot::color> octahedron_compound;

        ot::cell_set full_intersection;
        ot::cell_set octahedron_intersection;
    };

    const std::array<std::array<ot::octa_edge, 3>, triangle_count>& triangle_corner_table()
    {
        static constexpr std::array<std::array<ot::octa_edge, 3>, triangle_count> values = { {
            {
                ot::octa_edge::north_west_top,
                ot::octa_edge::east,
                ot::octa_edge::south_west_bottom
            },
            {
                ot::octa_edge::south_west_top,
                ot::octa_edge::south_east_bottom,
                ot::octa_edge::north
            },
            {
                ot::octa_edge::north_east_top,
                ot::octa_edge::south,
                ot::octa_edge::north_west_bottom
            },
            {
                ot::octa_edge::south_east_top,
                ot::octa_edge::west,
                ot::octa_edge::north_east_bottom
            }
        } };

        return values;
    }

    const std::array<ot::color, triangle_count>& triangle_color_table()
    {
        static const std::array<ot::color, triangle_count> values = {
            ot::color(0x1565c0),
            ot::color(0x2e7d32),
            ot::color(0xf9a825),
            ot::color(0xc62828)
        };

        return values;
    }

    ot::cell scaled_edge_cell(ot::octa_edge edge, int scale)
    {
        const ot::cell offset = ot::get_edge_neighbor(ot::cell(0, 0, 0), edge);
        return ot::cell(
            scale * offset.column,
            scale * offset.row,
            scale * offset.layer);
    }

    bool point_in_triangle(
        const ot::vector3& p,
        const ot::vector3& a,
        const ot::vector3& b,
        const ot::vector3& c,
        double eps = 1e-8)
    {
        const ot::vector3 n = (b - a).cross(c - a);
        const double n2 = n.squaredNorm();

        if (n2 < eps) {
            return false;
        }

        const ot::vector3 unit_n = n / std::sqrt(n2);
        if (std::abs((p - a).dot(unit_n)) > eps) {
            return false;
        }

        const ot::vector3 v0 = b - a;
        const ot::vector3 v1 = c - a;
        const ot::vector3 v2 = p - a;

        const double d00 = v0.dot(v0);
        const double d01 = v0.dot(v1);
        const double d11 = v1.dot(v1);
        const double d20 = v2.dot(v0);
        const double d21 = v2.dot(v1);

        const double denom = d00 * d11 - d01 * d01;
        if (std::abs(denom) < eps) {
            return false;
        }

        const double u = (d11 * d20 - d01 * d21) / denom;
        const double v = (d00 * d21 - d01 * d20) / denom;

        return u >= -eps && v >= -eps && u + v <= 1.0 + eps;
    }

    ot::cell_set make_cuboct_triangle_octas(
        const std::array<ot::octa_edge, 3>& corners,
        int scale)
    {
        const std::array<ot::cell, 3> corner_cells = {
            scaled_edge_cell(corners[0], scale),
            scaled_edge_cell(corners[1], scale),
            scaled_edge_cell(corners[2], scale)
        };

        const ot::vector3 a = corner_cells[0].center();
        const ot::vector3 b = corner_cells[1].center();
        const ot::vector3 c = corner_cells[2].center();

        int min_column = corner_cells[0].column;
        int max_column = corner_cells[0].column;
        int min_row = corner_cells[0].row;
        int max_row = corner_cells[0].row;
        int min_layer = corner_cells[0].layer;
        int max_layer = corner_cells[0].layer;

        for (const ot::cell& cell : corner_cells) {
            min_column = std::min(min_column, cell.column);
            max_column = std::max(max_column, cell.column);
            min_row = std::min(min_row, cell.row);
            max_row = std::max(max_row, cell.row);
            min_layer = std::min(min_layer, cell.layer);
            max_layer = std::max(max_layer, cell.layer);
        }

        ot::cell_set output;

        for (int column = min_column - 1; column <= max_column + 1; ++column) {
            for (int row = min_row - 1; row <= max_row + 1; ++row) {
                for (int layer = min_layer - 1; layer <= max_layer + 1; ++layer) {
                    const ot::cell cell(column, row, layer);
                    if (point_in_triangle(cell.center(), a, b, c)) {
                        output.add_cell(cell);
                    }
                }
            }
        }

        return output;
    }

    ot::cell_set get_pairwise_intersection(
        const std::array<ot::cell_set, triangle_count>& triangles,
        std::size_t count = triangle_count)
    {
        if (count > triangles.size()) {
            throw std::invalid_argument("Invalid triangle count");
        }

        ot::cell_set output;

        for (std::size_t i = 0; i < count; ++i) {
            for (std::size_t j = i + 1; j < count; ++j) {
                output.add_cells(triangles[i].intersect(triangles[j]));
            }
        }

        return output;
    }

    ot::colored_cell_set<ot::color> make_colored_compound(
        const std::array<ot::cell_set, triangle_count>& triangles,
        std::size_t count = triangle_count)
    {
        if (count > triangles.size()) {
            throw std::invalid_argument("Invalid triangle count");
        }

        ot::colored_cell_set<ot::color> output;
        const auto& colors = triangle_color_table();

        for (std::size_t i = 0; i < count; ++i) {
            output.add_cells(triangles[i].cells(), colors[i]);
        }

        return output;
    }

    ot::colored_cell_set<ot::color> color_cells_by_intersection(
        const ot::cell_set& cells,
        const ot::cell_set& intersection)
    {
        constexpr std::uint32_t white = 0xffffff;
        constexpr std::uint32_t red = 0xff0000;

        ot::colored_cell_set<ot::color> output;
        for (const ot::cell& cell : cells.cells()) {
            output.add_cell(
                cell,
                intersection.has_cell(cell) ? ot::color(red) : ot::color(white));
        }

        return output;
    }

    bool cell_is_in_intersection(
        const ot::cell& cell,
        const ot::cell_set& intersection)
    {
        if (intersection.has_cell(cell)) {
            return true;
        }

        // A split octahedron may now be represented by a pyramid, while
        // the intersection set still contains the original octahedron.
        return cell.type == ot::cell_type::pyramid &&
            intersection.has_cell(cell.to_octahedron());
    }

    ot::colored_cell_set<ot::color> color_cells_by_intersection_with_pyramids(
        const ot::cell_set& cells,
        const ot::cell_set& intersection)
    {
        constexpr std::uint32_t white = 0xffffff;
        constexpr std::uint32_t red = 0xff0000;

        ot::colored_cell_set<ot::color> output;
        for (const ot::cell& cell : cells.cells()) {
            output.add_cell(
                cell,
                cell_is_in_intersection(cell, intersection)
                ? ot::color(red)
                : ot::color(white));
        }

        return output;
    }

    cuboct_triangle_model make_cuboct_triangle_model(int scale)
    {
        cuboct_triangle_model output;
        const auto& triangle_corners = triangle_corner_table();

        for (std::size_t i = 0; i < triangle_corners.size(); ++i) {
            ot::cell_set triangle_octas =
                make_cuboct_triangle_octas(triangle_corners[i], scale);

            triangle_octas.remove_cell(global_center_octahedron());

            output.octahedron_triangles[i] = triangle_octas;

            output.full_triangles[i] =
                ot::fill_tetrahedra( output.octahedron_triangles[i]);
        }

        output.full_compound =
            make_colored_compound(output.full_triangles);

        output.octahedron_compound =
            make_colored_compound(output.octahedron_triangles);

        output.full_intersection =
            get_pairwise_intersection(output.full_triangles);

        output.octahedron_intersection =
            get_pairwise_intersection(output.octahedron_triangles);

        return output;
    }

    std::string output_path(
        const std::filesystem::path& output_dir,
        const std::string& filename)
    {
        return (output_dir / filename).string();
    }

    void write_cell_set(
        const ot::cell_set& cells,
        const std::filesystem::path& output_dir,
        const std::string& filename)
    {
        ot::cell_set_polyhedron poly(cells);
        poly.write_to_ply_file(output_path(output_dir, filename));
        std::println("wrote {}", filename);
    }

    void write_colored_cell_set(
        const ot::colored_cell_set<ot::color>& cells,
        const std::filesystem::path& output_dir,
        const std::string& filename,
        float buffer_percent = 1.0f)
    {
        ot::colored_cell_set_polyhedron poly(cells, buffer_percent);
        poly.write_to_ply_file(output_path(output_dir, filename));
        std::println("wrote {}", filename);
    }

    struct triangle_split_geometry {
        std::array<ot::cell, 3> corner_octas;
        std::array<ot::vector3, 3> corner_centers;
        std::array<ot::vector3, 3> side_centers;
    };

    double squared_distance_to_segment(
        const ot::vector3& point,
        const ot::vector3& a,
        const ot::vector3& b)
    {
        const ot::vector3 ab = b - a;
        const double ab2 = ab.squaredNorm();

        if (ab2 <= 1e-12) {
            return (point - a).squaredNorm();
        }

        const double t = std::clamp((point - a).dot(ab) / ab2, 0.0, 1.0);
        const ot::vector3 closest = a + t * ab;
        return (point - closest).squaredNorm();
    }

    int count_octahedron_edge_neighbors_in_set(
        const ot::cell& octahedron,
        const ot::cell_set& octahedra)
    {
        if (octahedron.type != ot::cell_type::octahedron) {
            throw std::invalid_argument("Expected octahedron");
        }

        int count = 0;
        for (const ot::cell& neighbor : octahedron.edge_neighbors()) {
            if (neighbor.type == ot::cell_type::octahedron &&
                octahedra.has_cell(neighbor)) {
                ++count;
            }
        }

        return count;
    }

    std::array<ot::cell, 3> get_triangle_corner_octas(const ot::cell_set& triangle)
    {
        const ot::cell_set octahedra(triangle.octahedra());

        int min_neighbors = std::numeric_limits<int>::max();
        std::vector<ot::cell> candidates;

        for (const ot::cell& octahedron : octahedra.cells()) {
            const int count =
                count_octahedron_edge_neighbors_in_set(octahedron, octahedra);

            if (count < min_neighbors) {
                min_neighbors = count;
                candidates.clear();
                candidates.push_back(octahedron);
            }
            else if (count == min_neighbors) {
                candidates.push_back(octahedron);
            }
        }

        if (candidates.size() != 3) {
            throw std::runtime_error(
                "Expected exactly three triangle-corner octahedra");
        }

        return { candidates[0], candidates[1], candidates[2] };
    }

    triangle_split_geometry get_triangle_split_geometry(const ot::cell_set& triangle)
    {
        triangle_split_geometry output;
        output.corner_octas = get_triangle_corner_octas(triangle);

        for (std::size_t i = 0; i < output.corner_octas.size(); ++i) {
            output.corner_centers[i] = output.corner_octas[i].center();
        }

        // Piece 0 corresponds to side between corners 1 and 2.
        // Piece 1 corresponds to side between corners 2 and 0.
        // Piece 2 corresponds to side between corners 0 and 1.
        output.side_centers[0] =
            0.5 * (output.corner_centers[1] + output.corner_centers[2]);

        output.side_centers[1] =
            0.5 * (output.corner_centers[2] + output.corner_centers[0]);

        output.side_centers[2] =
            0.5 * (output.corner_centers[0] + output.corner_centers[1]);

        return output;
    }

    std::vector<std::size_t> closest_indices(
        const ot::vector3& point,
        const std::array<ot::vector3, 3>& targets,
        double eps = 1e-8)
    {
        double best = std::numeric_limits<double>::max();
        std::vector<std::size_t> output;

        for (std::size_t i = 0; i < targets.size(); ++i) {
            const double d = (point - targets[i]).squaredNorm();

            if (d < best - eps) {
                best = d;
                output.clear();
                output.push_back(i);
            }
            else if (std::abs(d - best) <= eps) {
                output.push_back(i);
            }
        }

        return output;
    }

    std::size_t closest_side_to_tetrahedron(
        const ot::cell& tetrahedron,
        const triangle_split_geometry& geometry)
    {
        const ot::vector3 center = tetrahedron.center();

        const std::array<double, 3> side_distances = {
            squared_distance_to_segment(
                center,
                geometry.corner_centers[1],
                geometry.corner_centers[2]),

            squared_distance_to_segment(
                center,
                geometry.corner_centers[2],
                geometry.corner_centers[0]),

            squared_distance_to_segment(
                center,
                geometry.corner_centers[0],
                geometry.corner_centers[1])
        };

        double best = std::numeric_limits<double>::max();
        std::vector<std::size_t> closest;

        for (std::size_t i = 0; i < side_distances.size(); ++i) {
            const double d = side_distances[i];

            if (d < best - 1e-8) {
                best = d;
                closest.clear();
                closest.push_back(i);
            }
            else if (std::abs(d - best) <= 1e-8) {
                closest.push_back(i);
            }
        }

        if (closest.size() == 1) {
            return closest.front();
        }

        // Deterministic fallback for tetrahedra exactly on a side bisector:
        // choose among the tied sides by nearest side-center.
        double best_center_distance = std::numeric_limits<double>::max();
        std::size_t best_index = closest.front();

        for (std::size_t index : closest) {
            const double d = (center - geometry.side_centers[index]).squaredNorm();
            if (d < best_center_distance) {
                best_center_distance = d;
                best_index = index;
            }
        }

        return best_index;
    }

    ot::vector3 pyramid_geometric_center(const ot::cell& pyramid)
    {
        if (pyramid.type != ot::cell_type::pyramid) {
            throw std::invalid_argument("Expected pyramid");
        }

        ot::vector3 sum = ot::vector3::Zero();

        for (ot::octa_vert vert : ot::pyramid_vertices(pyramid.pyramid_tip())) {
            sum += ot::get_vertex_location(
                ot::get_octahedron_vertex(pyramid.to_octahedron(), vert));
        }

        return sum / 5.0;
    }

    const std::array<std::array<ot::octa_vert, 2>, 3>& opposite_pyramid_splits()
    {
        static constexpr std::array<std::array<ot::octa_vert, 2>, 3> values = { {
            { ot::octa_vert::top, ot::octa_vert::bottom },
            { ot::octa_vert::north_east, ot::octa_vert::south_west },
            { ot::octa_vert::north_west, ot::octa_vert::south_east }
        } };

        return values;
    }

    void split_octahedron_between_two_side_pieces(
        const ot::cell& octahedron,
        std::size_t side_a,
        std::size_t side_b,
        const triangle_split_geometry& geometry,
        std::array<ot::cell_set, 3>& pieces)
    {
        if (octahedron.type != ot::cell_type::octahedron) {
            throw std::invalid_argument("Expected octahedron");
        }

        struct candidate_split {
            ot::cell pyramid_for_a;
            ot::cell pyramid_for_b;
            double score = std::numeric_limits<double>::max();
        };

        std::optional<candidate_split> best;

        for (const auto& split : opposite_pyramid_splits()) {
            const ot::cell p0 =
                octahedron.to_pyramid(split[0]);

            const ot::cell p1 =
                octahedron.to_pyramid(split[1]);

            const ot::vector3 c0 = pyramid_geometric_center(p0);
            const ot::vector3 c1 = pyramid_geometric_center(p1);

            const double p0_to_a =
                (c0 - geometry.side_centers[side_a]).squaredNorm();

            const double p1_to_a =
                (c1 - geometry.side_centers[side_a]).squaredNorm();

            const bool p0_goes_to_a = p0_to_a <= p1_to_a;

            const ot::cell pyramid_for_a = p0_goes_to_a ? p0 : p1;
            const ot::cell pyramid_for_b = p0_goes_to_a ? p1 : p0;

            const double score =
                (pyramid_geometric_center(pyramid_for_a) -
                    geometry.side_centers[side_a]).squaredNorm() +
                (pyramid_geometric_center(pyramid_for_b) -
                    geometry.side_centers[side_b]).squaredNorm();

            if (!best || score < best->score) {
                best = candidate_split{
                    pyramid_for_a,
                    pyramid_for_b,
                    score
                };
            }
        }

        if (!best) {
            throw std::runtime_error("Could not split octahedron into pyramids");
        }

        pieces[side_a].add_cell(best->pyramid_for_a);
        pieces[side_b].add_cell(best->pyramid_for_b);
    }

    std::array<ot::cell_set, 3> split_triangle_into_three_pieces(
        const ot::cell_set& triangle)
    {
        const triangle_split_geometry geometry =
            get_triangle_split_geometry(triangle);

        std::array<ot::cell_set, 3> pieces;

        for (const ot::cell& cell : triangle.cells()) {
            if (cell.type == ot::cell_type::tetrahedron) {
                pieces[closest_side_to_tetrahedron(cell, geometry)].add_cell(cell);
                continue;
            }

            if (cell.type == ot::cell_type::pyramid) {
                const std::vector<std::size_t> closest =
                    closest_indices(pyramid_geometric_center(cell), geometry.side_centers);

                pieces[closest.front()].add_cell(cell);
                continue;
            }

            if (cell.type != ot::cell_type::octahedron) {
                throw std::runtime_error("Unknown cell type");
            }

            // This is the deliberate central void.
            if (cell == global_center_octahedron()) {
                continue;
            }

            const std::vector<std::size_t> closest =
                closest_indices(cell.center(), geometry.side_centers);

            if (closest.size() == 1) {
                pieces[closest.front()].add_cell(cell);
                continue;
            }

            if (closest.size() == 2) {
                split_octahedron_between_two_side_pieces(
                    cell,
                    closest[0],
                    closest[1],
                    geometry,
                    pieces);
                continue;
            }

            throw std::runtime_error(
                "An octahedron was equally closest to all three side centers");
        }

        return pieces;
    }

    ot::colored_cell_set<ot::color> color_triangle_pieces(
        const std::array<ot::cell_set, 3>& pieces)
    {
        const std::array<ot::color, 3> colors = {
            ot::color(0xff1744),
            ot::color(0x00c853),
            ot::color(0x2962ff)
        };

        ot::colored_cell_set<ot::color> output;

        for (std::size_t i = 0; i < pieces.size(); ++i) {
            output.add_cells(pieces[i].cells(), colors[i]);
        }

        return output;
    }
}

int main()
{
    const std::filesystem::path output_dir = "C:\\test\\octet2";
    std::filesystem::create_directories(output_dir);

    const cuboct_triangle_model model =
        make_cuboct_triangle_model(4);

    // 1. Full compound: octahedra + filled tetrahedra.
    write_colored_cell_set(
        model.full_compound,
        output_dir,
        "cuboct_triangle_compound.ply");

    // 2. Same compound, but octahedra only.
    write_colored_cell_set(
        model.octahedron_compound,
        output_dir,
        "cuboct_triangle_compound_octas.ply");

    // 3. Pairwise intersections among the four full triangles.
    write_cell_set(
        model.full_intersection,
        output_dir,
        "intersections.ply");

    // 4. One full triangle: intersection cells red, non-intersection cells white.
    write_colored_cell_set(
        color_cells_by_intersection(
            model.full_triangles[0],
            model.full_intersection),
        output_dir,
        "single_triangle_intersections_colored.ply",
        diagnostic_buffer_percent);

    // 5. Same as 4, but octahedra only.
    write_colored_cell_set(
        color_cells_by_intersection(
            model.octahedron_triangles[0],
            model.octahedron_intersection),
        output_dir,
        "single_triangle_intersections_colored_octas.ply",
        diagnostic_buffer_percent);

    // 6. Same as the full compound, but only the first two triangles.
    write_colored_cell_set(
        make_colored_compound(model.full_triangles, 2),
        output_dir,
        "cuboct_triangle_compound_two_triangles.ply");

    const std::array<ot::cell_set, 3> triangle_pieces =
        split_triangle_into_three_pieces(model.full_triangles[0]);

    // 7. One triangle split into its three colored pieces.
    write_colored_cell_set(
        color_triangle_pieces(triangle_pieces),
        output_dir,
        "single_triangle_three_piece_split.ply",
        diagnostic_buffer_percent);

    // 8. First piece from the three-piece split.
    write_cell_set(
        triangle_pieces[0],
        output_dir,
        "single_triangle_piece_0.ply");

    // 9. Same first piece, but with tetrahedra removed.
    write_cell_set(
        ot::cell_set(triangle_pieces[0].octahedra_and_pyramids()),
        output_dir,
        "single_triangle_piece_0_octas_and_pyramids.ply");

    // 10. Same first piece, but with intersection cells red and non-intersection cells white.
    write_colored_cell_set(
        color_cells_by_intersection_with_pyramids(
            triangle_pieces[0],
            model.full_intersection),
        output_dir,
        "single_triangle_piece_0_intersections_colored.ply",
        diagnostic_buffer_percent);

    write_cell_set(
        ot::create_piece(4, 1),
        output_dir,
        "piece_work.ply");
}