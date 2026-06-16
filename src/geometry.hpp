#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ot {

    using vector2 = Eigen::Vector2d;
    using vector3 = Eigen::Vector3d;
    using vector4 = Eigen::Vector4d;
    using matrix3 = Eigen::Matrix3d;
    using matrix4 = Eigen::Matrix4d;

    inline constexpr double default_epsilon = 1e-12;
    inline constexpr double ray_triangle_epsilon = 1e-8;

    namespace detail {

    inline bool nearly_zero(double value, double epsilon = default_epsilon)
    {
        return std::abs(value) <= epsilon;
    }

    inline vector3 normalized_or_zero(const vector3& value, double epsilon = default_epsilon)
    {
        const double length = value.norm();
        if (length <= epsilon) {
            return vector3::Zero();
        }

        return value / length;
    }

    inline int dominant_axis(const vector3& value)
    {
        const vector3 abs_value = value.cwiseAbs();
        if (abs_value.x() >= abs_value.y() && abs_value.x() >= abs_value.z()) {
            return 0;
        }

        if (abs_value.y() >= abs_value.z()) {
            return 1;
        }

        return 2;
    }

    inline int projection_axis_to_drop(const vector3& normal)
    {
        return dominant_axis(normal);
    }

    inline vector2 project_to_2d(const vector3& value, int axis_to_drop)
    {
        switch (axis_to_drop) {
        case 0:
            return { value.y(), value.z() };
        case 1:
            return { value.x(), value.z() };
        default:
            return { value.x(), value.y() };
        }
    }

    inline double orient_2d(const vector2& a, const vector2& b, const vector2& c)
    {
        const vector2 ab = b - a;
        const vector2 ac = c - a;
        return ab.x() * ac.y() - ab.y() * ac.x();
    }

    inline bool point_on_segment_2d(
        const vector2& point,
        const vector2& a,
        const vector2& b,
        double epsilon)
    {
        if (std::abs(orient_2d(a, b, point)) > epsilon) {
            return false;
        }

        return point.x() >= std::min(a.x(), b.x()) - epsilon &&
               point.x() <= std::max(a.x(), b.x()) + epsilon &&
               point.y() >= std::min(a.y(), b.y()) - epsilon &&
               point.y() <= std::max(a.y(), b.y()) + epsilon;
    }

    inline bool segments_intersect_2d(
        const vector2& a0,
        const vector2& a1,
        const vector2& b0,
        const vector2& b1,
        double epsilon)
    {
        const double o1 = orient_2d(a0, a1, b0);
        const double o2 = orient_2d(a0, a1, b1);
        const double o3 = orient_2d(b0, b1, a0);
        const double o4 = orient_2d(b0, b1, a1);

        if (((o1 > epsilon && o2 < -epsilon) || (o1 < -epsilon && o2 > epsilon)) &&
            ((o3 > epsilon && o4 < -epsilon) || (o3 < -epsilon && o4 > epsilon))) {
            return true;
        }

        return point_on_segment_2d(b0, a0, a1, epsilon) ||
               point_on_segment_2d(b1, a0, a1, epsilon) ||
               point_on_segment_2d(a0, b0, b1, epsilon) ||
               point_on_segment_2d(a1, b0, b1, epsilon);
    }

    inline bool point_in_triangle_2d(
        const vector2& point,
        const vector2& a,
        const vector2& b,
        const vector2& c,
        double epsilon)
    {
        const double o0 = orient_2d(a, b, point);
        const double o1 = orient_2d(b, c, point);
        const double o2 = orient_2d(c, a, point);

        const bool has_negative = o0 < -epsilon || o1 < -epsilon || o2 < -epsilon;
        const bool has_positive = o0 > epsilon || o1 > epsilon || o2 > epsilon;

        return !(has_negative && has_positive);
    }

    inline bool coplanar_tri_tri(
        const vector3& normal,
        const vector3& v0,
        const vector3& v1,
        const vector3& v2,
        const vector3& u0,
        const vector3& u1,
        const vector3& u2,
        double epsilon)
    {
        const int axis_to_drop = projection_axis_to_drop(normal);

        const std::array<vector2, 3> v = {
            project_to_2d(v0, axis_to_drop),
            project_to_2d(v1, axis_to_drop),
            project_to_2d(v2, axis_to_drop)
        };

        const std::array<vector2, 3> u = {
            project_to_2d(u0, axis_to_drop),
            project_to_2d(u1, axis_to_drop),
            project_to_2d(u2, axis_to_drop)
        };

        for (int i = 0; i < 3; ++i) {
            const int inext = (i + 1) % 3;
            for (int j = 0; j < 3; ++j) {
                const int jnext = (j + 1) % 3;
                if (segments_intersect_2d(v[i], v[inext], u[j], u[jnext], epsilon)) {
                    return true;
                }
            }
        }

        return point_in_triangle_2d(v[0], u[0], u[1], u[2], epsilon) ||
               point_in_triangle_2d(u[0], v[0], v[1], v[2], epsilon);
    }

    inline double cleaned_signed_distance(double value, double epsilon)
    {
        return std::abs(value) < epsilon ? 0.0 : value;
    }

    inline bool all_same_strict_side(double d0, double d1, double d2)
    {
        return (d0 > 0.0 && d1 > 0.0 && d2 > 0.0) ||
               (d0 < 0.0 && d1 < 0.0 && d2 < 0.0);
    }

    inline std::optional<std::pair<double, double>> triangle_plane_line_interval(
        const std::array<vector3, 3>& triangle,
        const std::array<double, 3>& signed_distances,
        int projection_axis,
        double epsilon)
    {
        std::vector<double> values;
        values.reserve(4);

        auto add_value = [&](double value) {
            for (double existing : values) {
                if (std::abs(existing - value) <= epsilon) {
                    return;
                }
            }
            values.push_back(value);
        };

        for (int i = 0; i < 3; ++i) {
            if (signed_distances[i] == 0.0) {
                add_value(triangle[i][projection_axis]);
            }
        }

        for (int i = 0; i < 3; ++i) {
            const int j = (i + 1) % 3;
            const double di = signed_distances[i];
            const double dj = signed_distances[j];

            if (di * dj < 0.0) {
                const double t = di / (di - dj);
                const vector3 point = triangle[i] + t * (triangle[j] - triangle[i]);
                add_value(point[projection_axis]);
            }
        }

        if (values.empty()) {
            return std::nullopt;
        }

        const auto [min_it, max_it] = std::minmax_element(values.begin(), values.end());
        return std::pair<double, double>(*min_it, *max_it);
    }

    } // namespace detail

    struct plane {
        vector3 normal = vector3::UnitZ();
        double d = 0.0;

        plane() = default;

        plane(const vector3& normal_, double d_)
            : normal(normal_)
            , d(d_)
        {
        }

        explicit plane(const vector4& value)
            : normal(value.x(), value.y(), value.z())
            , d(value.w())
        {
        }

        plane(double a, double b, double c, double d_)
            : normal(a, b, c)
            , d(d_)
        {
        }

        plane(const vector3& a, const vector3& b, const vector3& c)
        {
            normal = detail::normalized_or_zero((b - a).cross(c - a));
            d = -normal.dot(a);
        }

        [[nodiscard]] vector4 coefficients() const
        {
            return { normal.x(), normal.y(), normal.z(), d };
        }

        [[nodiscard]] double dot(const vector4& value) const
        {
            return coefficients().dot(value);
        }

        [[nodiscard]] double dot_coordinate(const vector3& value) const
        {
            return normal.dot(value) + d;
        }

        [[nodiscard]] double dot_normal(const vector3& value) const
        {
            return normal.dot(value);
        }

        void normalize(double epsilon = default_epsilon)
        {
            const double length = normal.norm();
            if (length <= epsilon) {
                normal = vector3::Zero();
                d = 0.0;
                return;
            }

            normal /= length;
            d /= length;
        }

        [[nodiscard]] static plane normalized(plane value, double epsilon = default_epsilon)
        {
            value.normalize(epsilon);
            return value;
        }

        [[nodiscard]] static plane transform(const plane& value, const matrix4& transform)
        {
            const vector4 transformed = transform.inverse().transpose() * value.coefficients();
            return plane(transformed);
        }
    };

    inline bool operator==(const plane& lhs, const plane& rhs)
    {
        return lhs.normal == rhs.normal && lhs.d == rhs.d;
    }

    inline bool operator!=(const plane& lhs, const plane& rhs)
    {
        return !(lhs == rhs);
    }

    inline double classify_point(const vector3& point, const plane& value)
    {
        return value.dot_coordinate(point);
    }

    inline double perpendicular_distance(const vector3& point, const plane& value, double epsilon = default_epsilon)
    {
        const double denominator = value.normal.norm();
        if (denominator <= epsilon) {
            return 0.0;
        }

        return std::abs(value.normal.dot(point) + value.d) / denominator;
    }

    inline matrix4 rotate_around_line(const vector3& p1, const vector3& p2, double angle_radians)
    {
        const vector3 axis = detail::normalized_or_zero(p2 - p1);
        if (axis.isZero(default_epsilon)) {
            return matrix4::Identity();
        }

        const Eigen::Translation3d translate_to_origin(-p1);
        const Eigen::AngleAxisd rotate_around_axis(angle_radians, axis);
        const Eigen::Translation3d translate_back(p1);

        return (translate_back * rotate_around_axis * translate_to_origin).matrix();
    }

    inline vector3 transform_point(const vector3& point, const matrix4& transform, double epsilon = default_epsilon)
    {
        const vector4 transformed = transform * vector4(point.x(), point.y(), point.z(), 1.0);
        if (std::abs(transformed.w()) > epsilon && std::abs(transformed.w() - 1.0) > epsilon) {
            return transformed.head<3>() / transformed.w();
        }

        return transformed.head<3>();
    }

    inline std::vector<vector3> apply(const matrix4& transform, const std::vector<vector3>& values)
    {
        std::vector<vector3> transformed;
        transformed.reserve(values.size());

        for (const vector3& value : values) {
            transformed.push_back(transform_point(value, transform));
        }

        return transformed;
    }

    inline vector3 normalize(const vector3& value, double epsilon = default_epsilon)
    {
        return detail::normalized_or_zero(value, epsilon);
    }

    inline vector3 transform_vector(const vector3& value, const matrix4& transform)
    {
        return transform.block<3, 3>(0, 0) * value;
    }

    inline vector3 transform_normal(
        const vector3& normal,
        const matrix4& transform,
        bool normalize_result = true,
        double epsilon = default_epsilon)
    {
        const vector3 transformed = transform.block<3, 3>(0, 0).inverse().transpose() * normal;

        if (!normalize_result) {
            return transformed;
        }

        return detail::normalized_or_zero(transformed, epsilon);
    }

    inline vector3 get_face_normal(const std::vector<vector3>& face, double epsilon = default_epsilon)
    {
        if (face.size() < 3) {
            throw std::invalid_argument("get_face_normal requires at least three points");
        }

        return detail::normalized_or_zero((face[1] - face[0]).cross(face[2] - face[0]), epsilon);
    }

    inline vector3 get_face_centroid(const std::vector<vector3>& face)
    {
        if (face.empty()) {
            throw std::invalid_argument("get_face_centroid requires at least one point");
        }

        vector3 sum = vector3::Zero();
        for (const vector3& point : face) {
            sum += point;
        }

        return sum / static_cast<double>(face.size());
    }

    inline vector3 lerp(const vector3& p1, const vector3& p2, double t)
    {
        return p1 + (p2 - p1) * t;
    }

    inline bool line_line_intersect(
        const vector3& p1,
        const vector3& p2,
        const vector3& p3,
        const vector3& p4,
        double epsilon,
        vector3& pa,
        vector3& pb,
        double& mua,
        double& mub)
    {
        const vector3 p13 = p1 - p3;
        const vector3 p43 = p4 - p3;
        if (std::abs(p43.x()) < epsilon && std::abs(p43.y()) < epsilon && std::abs(p43.z()) < epsilon) {
            return false;
        }

        const vector3 p21 = p2 - p1;
        if (std::abs(p21.x()) < epsilon && std::abs(p21.y()) < epsilon && std::abs(p21.z()) < epsilon) {
            return false;
        }

        const double d1343 = p13.dot(p43);
        const double d4321 = p43.dot(p21);
        const double d1321 = p13.dot(p21);
        const double d4343 = p43.dot(p43);
        const double d2121 = p21.dot(p21);

        const double denom = d2121 * d4343 - d4321 * d4321;
        if (std::abs(denom) < epsilon) {
            return false;
        }

        const double numer = d1343 * d4321 - d1321 * d4343;
        mua = numer / denom;
        mub = (d1343 + d4321 * mua) / d4343;

        pa = p1 + mua * p21;
        pb = p3 + mub * p43;
        return true;
    }

    inline std::optional<vector3> line_line_intersect(
        const vector3& p1,
        const vector3& p2,
        const vector3& p3,
        const vector3& p4,
        double epsilon = default_epsilon)
    {
        vector3 pa = vector3::Zero();
        vector3 pb = vector3::Zero();
        double mua = 0.0;
        double mub = 0.0;

        if (!line_line_intersect(p1, p2, p3, p4, epsilon, pa, pb, mua, mub)) {
            return std::nullopt;
        }

        return (pa + pb) / 2.0;
    }

    inline std::optional<double> ray_triangle_intersection(
        const vector3& origin,
        const vector3& direction,
        const vector3& v0,
        const vector3& v1,
        const vector3& v2,
        double epsilon = ray_triangle_epsilon,
        double t_min = 0.0)
    {
        const vector3 e1 = v1 - v0;
        const vector3 e2 = v2 - v0;

        const vector3 pvec = direction.cross(e2);
        const double det = e1.dot(pvec);

        if (std::abs(det) < epsilon) {
            return std::nullopt;
        }

        const double inv_det = 1.0 / det;
        const vector3 tvec = origin - v0;
        const double u = tvec.dot(pvec) * inv_det;

        if (u < 0.0 || u > 1.0) {
            return std::nullopt;
        }

        const vector3 qvec = tvec.cross(e1);
        const double v = direction.dot(qvec) * inv_det;

        if (v < 0.0 || u + v > 1.0) {
            return std::nullopt;
        }

        const double t = e2.dot(qvec) * inv_det;
        if (t < t_min) {
            return std::nullopt;
        }

        return t;
    }

    inline double ray_triangle_intersection_or_zero(
        const vector3& origin,
        const vector3& direction,
        const vector3& v0,
        const vector3& v1,
        const vector3& v2,
        double epsilon = ray_triangle_epsilon,
        double t_min = 0.0)
    {
        return ray_triangle_intersection(origin, direction, v0, v1, v2, epsilon, t_min).value_or(0.0);
    }

    inline bool has_triangle_triangle_intersection(
        const vector3& v0,
        const vector3& v1,
        const vector3& v2,
        const vector3& u0,
        const vector3& u1,
        const vector3& u2,
        double epsilon = default_epsilon)
    {
        const vector3 e1 = v1 - v0;
        const vector3 e2 = v2 - v0;
        const vector3 n1 = e1.cross(e2);

        const vector3 f1 = u1 - u0;
        const vector3 f2 = u2 - u0;
        const vector3 n2 = f1.cross(f2);

        if (n1.squaredNorm() <= epsilon * epsilon || n2.squaredNorm() <= epsilon * epsilon) {
            return false;
        }

        const double d1 = -n1.dot(v0);
        const double du0 = detail::cleaned_signed_distance(n1.dot(u0) + d1, epsilon);
        const double du1 = detail::cleaned_signed_distance(n1.dot(u1) + d1, epsilon);
        const double du2 = detail::cleaned_signed_distance(n1.dot(u2) + d1, epsilon);

        if (detail::all_same_strict_side(du0, du1, du2)) {
            return false;
        }

        const double d2 = -n2.dot(u0);
        const double dv0 = detail::cleaned_signed_distance(n2.dot(v0) + d2, epsilon);
        const double dv1 = detail::cleaned_signed_distance(n2.dot(v1) + d2, epsilon);
        const double dv2 = detail::cleaned_signed_distance(n2.dot(v2) + d2, epsilon);

        if (detail::all_same_strict_side(dv0, dv1, dv2)) {
            return false;
        }

        const vector3 intersection_direction = n1.cross(n2);
        if (intersection_direction.squaredNorm() <= epsilon * epsilon) {
            if (du0 == 0.0 && du1 == 0.0 && du2 == 0.0 &&
                dv0 == 0.0 && dv1 == 0.0 && dv2 == 0.0) {
                return detail::coplanar_tri_tri(n1, v0, v1, v2, u0, u1, u2, epsilon);
            }

            return false;
        }

        const int projection_axis = detail::dominant_axis(intersection_direction);

        const std::array<vector3, 3> v_triangle = { v0, v1, v2 };
        const std::array<vector3, 3> u_triangle = { u0, u1, u2 };
        const std::array<double, 3> v_distances = { dv0, dv1, dv2 };
        const std::array<double, 3> u_distances = { du0, du1, du2 };

        const auto v_interval = detail::triangle_plane_line_interval(
            v_triangle,
            v_distances,
            projection_axis,
            epsilon);

        const auto u_interval = detail::triangle_plane_line_interval(
            u_triangle,
            u_distances,
            projection_axis,
            epsilon);

        if (!v_interval || !u_interval) {
            return false;
        }

        return v_interval->first <= u_interval->second + epsilon &&
               u_interval->first <= v_interval->second + epsilon;
    }

    inline bool are_coplanar(
        const vector3& v1,
        const vector3& v2,
        const vector3& v3,
        const vector3& v4,
        double epsilon = default_epsilon)
    {
        const double value = (v2 - v1).cross(v4 - v1).dot(v3 - v1);
        return std::abs(value) < epsilon;
    }

} // namespace geom3d
