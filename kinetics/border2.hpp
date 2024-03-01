/**
 * @file border2.hpp
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
#include <cstddef>
#include <cstdint>
#include <utility>

#include "common/dualism.hpp"
#include "common/utils/macros.hpp"
#include "math/curves/piecewise_curves/piecewise_linear_curve2.hpp"

namespace boyle {
namespace kinetics {

template <std::floating_point T>
class [[nodiscard]] Border2 {
  public:
    explicit Border2(
        std::uint64_t c_id, ::boyle::common::Chirality c_chirality,
        ::boyle::math::PiecewiseLinearCurve2<T> c_bound_line
    ) noexcept
        : id{c_id}, chirality{c_chirality}, bound_line{std::move(c_bound_line)} {}
    ENABLE_IMPLICIT_CONSTRUCTORS(Border2);
    virtual ~Border2() noexcept = 0;

    std::uint64_t id;
    ::boyle::common::Chirality chirality;
    ::boyle::math::PiecewiseLinearCurve2<T> bound_line;
};

template <std::floating_point T>
inline Border2<T>::~Border2() noexcept = default;

template <std::floating_point T>
class [[nodiscard]] HardBorder2 final : public Border2<T> {
  public:
    explicit HardBorder2(
        std::uint64_t c_id, ::boyle::common::Chirality c_chirality,
        ::boyle::math::PiecewiseLinearCurve2<T> c_bound_line
    ) noexcept
        : Border2<T>(c_id, c_chirality, std::move(c_bound_line)) {}
    ENABLE_IMPLICIT_CONSTRUCTORS(HardBorder2);
    ~HardBorder2() noexcept override = default;
};

template <std::floating_point T>
class [[nodiscard]] SoftBorder2 final : public Border2<T> {
  public:
    explicit SoftBorder2(
        std::uint64_t c_id, ::boyle::common::Chirality c_chirality,
        ::boyle::math::PiecewiseLinearCurve2<T> c_bound_line, T c_linear_weight,
        T c_quadratic_weight
    ) noexcept
        : Border2<T>(c_id, c_chirality, std::move(c_bound_line)), linear_weight{c_linear_weight},
          quadratic_weight{c_quadratic_weight} {}
    ENABLE_IMPLICIT_CONSTRUCTORS(SoftBorder2);
    ~SoftBorder2() noexcept override = default;

    T linear_weight;
    T quadratic_weight;
};

using HardBorder2f = HardBorder2<float>;
using HardBorder2d = HardBorder2<double>;

using SoftBorder2f = SoftBorder2<float>;
using SoftBorder2d = SoftBorder2<double>;

} // namespace kinetics
} // namespace boyle

namespace boost {
namespace serialization {

template <typename Archive, std::floating_point T>
[[using gnu: always_inline]]
inline void serialize(
    Archive& ar, ::boyle::kinetics::HardBorder2<T>& obj, const unsigned int version
) noexcept {
    ar& obj.id;
    ar& obj.chirality;
    ar& obj.bound_line;
    return;
}

template <typename Archive, std::floating_point T>
[[using gnu: always_inline]]
inline void serialize(
    Archive& ar, ::boyle::kinetics::SoftBorder2<T>& obj, const unsigned int version
) noexcept {
    ar& obj.id;
    ar& obj.chirality;
    ar& obj.bound_line;
    ar& obj.linear_weight;
    ar& obj.quadratic_weight;
    return;
}

} // namespace serialization
} // namespace boost
