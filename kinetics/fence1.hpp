/**
 * @file fence1.hpp
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
#include "math/functions/piecewise_functions/piecewise_linear_function1.hpp"

namespace boyle {
namespace kinetics {

template <std::floating_point T>
class [[nodiscard]] Fence1 {
  public:
    explicit Fence1(
        std::uint64_t c_id, ::boyle::common::Actio c_actio,
        ::boyle::math::PiecewiseLinearFunction1<T> c_bound_line
    ) noexcept
        : id{c_id}, actio{c_actio}, bound_line{std::move(c_bound_line)} {}
    ENABLE_IMPLICIT_CONSTRUCTORS(Fence1);
    virtual ~Fence1() noexcept = 0;

    std::uint64_t id;
    ::boyle::common::Actio actio;
    ::boyle::math::PiecewiseLinearFunction1<T> bound_line;
};

template <std::floating_point T>
inline Fence1<T>::~Fence1() noexcept = default;

template <std::floating_point T>
class [[nodiscard]] HardFence1 final : public Fence1<T> {
  public:
    explicit HardFence1(
        std::uint64_t c_id, ::boyle::common::Actio c_actio,
        ::boyle::math::PiecewiseLinearFunction1<T> c_bound_line
    ) noexcept
        : Fence1<T>(c_id, c_actio, std::move(c_bound_line)) {}
    ENABLE_IMPLICIT_CONSTRUCTORS(HardFence1);
    ~HardFence1() noexcept override = default;
};

template <std::floating_point T>
class [[nodiscard]] SoftFence1 final : public Fence1<T> {
  public:
    explicit SoftFence1(
        std::uint64_t c_id, ::boyle::common::Actio c_actio,
        ::boyle::math::PiecewiseLinearFunction1<T> c_bound_line, T c_linear_weight,
        T c_quadratic_weight
    ) noexcept
        : Fence1<T>(c_id, c_actio, std::move(c_bound_line)), linear_weight{c_linear_weight},
          quadratic_weight{c_quadratic_weight} {}
    ENABLE_IMPLICIT_CONSTRUCTORS(SoftFence1);
    ~SoftFence1() noexcept override = default;

    T linear_weight;
    T quadratic_weight;
};

using HardFence1f = HardFence1<float>;
using HardFence1d = HardFence1<double>;

using SoftFence1f = SoftFence1<float>;
using SoftFence1d = SoftFence1<double>;

} // namespace kinetics
} // namespace boyle

namespace boost {
namespace serialization {

template <typename Archive, std::floating_point T>
[[using gnu: always_inline]]
inline void serialize(
    Archive& ar, boyle::kinetics::HardFence1<T>& obj, const unsigned int version
) noexcept {
    ar& obj.id;
    ar& obj.actio;
    ar& obj.bound_line;
    return;
}

template <typename Archive, std::floating_point T>
[[using gnu: always_inline]]
inline void serialize(
    Archive& ar, boyle::kinetics::SoftFence1<T>& obj, const unsigned int version
) noexcept {
    ar& obj.id;
    ar& obj.actio;
    ar& obj.bound_line;
    ar& obj.linear_weight;
    ar& obj.quadratic_weight;
    return;
}

} // namespace serialization
} // namespace boost
