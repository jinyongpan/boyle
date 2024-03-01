/**
 * @file path2.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-06
 *
 * @copyright Copyright (c) 2023 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include <array>
#include <concepts>
#include <utility>
#include <vector>

#include "boost/serialization/access.hpp"

#include "common/utils/macros.hpp"
#include "math/curves/piecewise_curves/piecewise_quintic_curve2.hpp"
#include "math/duplet.hpp"
#include "math/vec2.hpp"

namespace boyle {
namespace kinetics {

template <std::floating_point T>
class [[nodiscard]] Path2 final {
    friend class boost::serialization::access;

  public:
    using BoundaryMode = typename ::boyle::math::PiecewiseQuinticCurve2<T>::BoundaryMode;

    [[using gnu: always_inline]] explicit Path2(
        std::vector<::boyle::math::Vec2<T>> anchor_points, T s0 = 0.0
    )
        : Path2(
              std::move(anchor_points),
              std::array<BoundaryMode, 2>{
                  BoundaryMode{2, ::boyle::math::Vec2<T>{0.0, 0.0}},
                  BoundaryMode{4, ::boyle::math::Vec2<T>{0.0, 0.0}}},
              std::array<BoundaryMode, 2>{
                  BoundaryMode{2, ::boyle::math::Vec2<T>{0.0, 0.0}},
                  BoundaryMode{4, ::boyle::math::Vec2<T>{0.0, 0.0}}},
              s0
          ) {}

    [[using gnu: always_inline]] explicit Path2(
        std::vector<::boyle::math::Vec2<T>> anchor_points, std::array<BoundaryMode, 2> b0,
        std::array<BoundaryMode, 2> bf, T s0 = 0.0
    )
        : curve_(anchor_points, b0, bf, s0) {}

    ENABLE_IMPLICIT_CONSTRUCTORS(Path2);

    ~Path2() noexcept = default;

    [[using gnu: pure, always_inline]]
    ::boyle::math::Vec2<T>
    operator()(T s) const noexcept {
        return curve_(s);
    }

    [[using gnu: pure, always_inline]]
    ::boyle::math::Vec2<T>
    operator()(T s, T l) const noexcept {
        return curve_(s, l);
    }

    [[using gnu: pure, always_inline]]
    ::boyle::math::Vec2<T>
    operator()(::boyle::math::SlDuplet<T> sl) const noexcept {
        return curve_(sl);
    }

    [[using gnu: pure, always_inline]]
    ::boyle::math::SlDuplet<T> inverse(::boyle::math::Vec2<T> point) const noexcept {
        return curve_.inverse(point);
    }

    [[using gnu: pure, always_inline]]
    ::boyle::math::SlDuplet<T> inverse(::boyle::math::Vec2<T> point, T start_s, T end_s)
        const noexcept {
        return curve_.inverse(point, start_s, end_s);
    }

    [[using gnu: pure, always_inline]]
    ::boyle::math::Vec2<T> tangent(T s) const noexcept {
        return curve_.tangent(s);
    }

    [[using gnu: pure, always_inline]]
    ::boyle::math::Vec2<T> orthnormal(T s) const noexcept {
        return curve_.orthonormal(s);
    }

    [[using gnu: pure, always_inline]]
    T curvature(T s) {
        return curve_.curvature(s);
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    std::vector<::boyle::math::Vec2<T>>
    operator()(const std::vector<T>& ss) const noexcept {
        return curve_(ss);
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    std::vector<::boyle::math::Vec2<T>> tangent(const std::vector<T>& ss) const noexcept {
        return curve_.tangent(ss);
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    std::vector<::boyle::math::Vec2<T>> orthonormal(const std::vector<T>& ss) const noexcept {
        return curve_.orthonormal(ss);
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    std::vector<::boyle::math::Vec2<T>>
    operator()(const std::vector<::boyle::math::SlDuplet<T>>& sls) const noexcept {
        return curve_(sls);
    }

  private:
    [[using gnu: always_inline]] explicit Path2(::boyle::math::PiecewiseQuinticCurve2<T> curve
    ) noexcept
        : curve_{std::move(curve)} {}

    template <typename Archive>
    [[using gnu: always_inline]]
    void serialize(Archive& ar) noexcept {
        ar& curve_;
        return;
    }

    ::boyle::math::PiecewiseQuinticCurve2<T> curve_;
};

using Path2f = Path2<float>;
using Path2d = Path2<double>;

} // namespace kinetics
} // namespace boyle
