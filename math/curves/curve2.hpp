/**
 * @file curve2.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-24
 *
 * @copyright Copyright (c) 2023 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>

#include "common/utils/macros.hpp"
#include "math/duplet.hpp"
#include "math/functions/function1.hpp"
#include "math/vec2.hpp"

namespace boyle {
namespace math {

template <typename DerivedCurve2, std::floating_point T>
class [[nodiscard]] Curve2 {
  public:
    using value_type = Vec2<T>;

    ENABLE_IMPLICIT_CONSTRUCTORS(Curve2);
    virtual ~Curve2() noexcept = 0;

    [[using gnu: pure, always_inline]]
    Vec2<T>
    operator()(T s) const noexcept {
        return underlying()->operator()(s);
    }

    [[using gnu: pure, always_inline]]
    Vec2<T> tangent(T s) const noexcept {
        return underlying()->tangent(s);
    }

    [[using gnu: pure, always_inline]]
    Vec2<T> inverse(Vec2<T> point) const noexcept {
        return underlying()->inverse(point);
    }

    [[using gnu: pure, always_inline]]
    Vec2<T> inverse(Vec2<T> point, T start_s, T end_s) const noexcept {
        return underlying()->inverse(point, start_s, end_s);
    }

    [[using gnu: pure, always_inline]]
    T minS() const noexcept {
        return underlying()->minS();
    }

    [[using gnu: pure, always_inline]]
    T maxS() const noexcept {
        return underlying()->maxS();
    }

    [[using gnu: pure, always_inline]]
    Vec2<T> orthonormal(T s) const noexcept {
        return tangent(s).selfRotateHalfPi();
    }

    [[using gnu: pure, always_inline]]
    Vec2<T>
    operator()(T s, T l) const noexcept {
        return operator()(s) + l * orthonormal(s);
    }

    [[using gnu: pure, always_inline]]
    Vec2<T>
    operator()(SlDuplet<T> sl) const noexcept {
        return operator()(sl.s) + sl.l * orthonormal(sl.s);
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    std::vector<Vec2<T>>
    operator()(const std::vector<T>& ss) const noexcept {
        const std::size_t size{ss.size()};
        std::vector<Vec2<T>> points;
        points.reserve(size);
        for (const T s : ss) {
            points.emplace_back(operator()(s));
        }
        return points;
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    std::vector<Vec2<T>> tangent(const std::vector<T>& ss) const noexcept {
        const std::size_t size{ss.size()};
        std::vector<Vec2<T>> tangents;
        tangents.reserve(size);
        for (const T s : ss) {
            tangents.emplace_back(tangent(s));
        }
        return tangents;
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    std::vector<Vec2<T>> orthonormal(const std::vector<T>& ss) const noexcept {
        const std::size_t size{ss.size()};
        std::vector<Vec2<T>> orthonormals;
        orthonormals.reserve(size);
        for (const T s : ss) {
            orthonormals.emplace_back(orthonormal(s));
        }
        return orthonormals;
    }

    [[using gnu: pure, always_inline]] [[nodiscard]]
    std::vector<Vec2<T>>
    operator()(const std::vector<SlDuplet<T>>& sls) const noexcept {
        const std::size_t size{sls.size()};
        std::vector<Vec2<T>> points;
        points.reserve(size);
        for (const T sl : sls) {
            points.emplace_back(operator()(sl));
        }
        return points;
    }

    [[using gnu: pure, always_inline]]
    Vec2<T> front() const noexcept {
        return operator()(minS());
    }

    [[using gnu: pure, always_inline]]
    Vec2<T> back() const noexcept {
        return operator()(maxS());
    }

  private:
    [[using gnu: pure, always_inline]]
    DerivedCurve2* underlying() noexcept {
        return static_cast<DerivedCurve2*>(this);
    }

    [[using gnu: pure, always_inline]]
    const DerivedCurve2* underlying() const noexcept {
        return static_cast<const DerivedCurve2*>(this);
    }
};

template <typename DerivedCurve2, std::floating_point T>
inline Curve2<DerivedCurve2, T>::~Curve2() noexcept = default;

} // namespace math
} // namespace boyle
