/**
 * @file vec3.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-08-02
 *
 * @copyright Copyright (c) 2023 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include <cmath>
#include <concepts>
#include <ostream>

#include "math/concepts.hpp"
#include "math/triplet.hpp"

namespace boyle {
namespace math {

template <std::floating_point T>
class Vec3 final : public Triplet<T> {
    using Triplet<T>::first_;
    using Triplet<T>::second_;
    using Triplet<T>::third_;

  public:
    using value_type = T;

    [[using gnu: always_inline]] Vec3() noexcept : Triplet<T>{} {}
    [[using gnu: always_inline]] constexpr Vec3(T cv) noexcept : Triplet<T>{cv} {}
    [[using gnu: always_inline]] constexpr Vec3(T cx, T cy, T cz) noexcept
        : Triplet<T>{cx, cy, cz} {}
    [[using gnu: always_inline]] constexpr Vec3(const Vec3& other) noexcept : Triplet<T>{other} {}
    [[using gnu: always_inline]]
    constexpr Vec3&
    operator=(const Vec3& other) noexcept {
        Triplet<T>::operator=(other);
        return *this;
    }
    ~Vec3() noexcept override {}

    template <std::floating_point U>
    [[using gnu: pure, always_inline]] constexpr operator Vec3<U>() noexcept {
        return Vec3<U>{static_cast<U>(x), static_cast<U>(y), static_cast<U>(z)};
    }

    [[using gnu: pure, always_inline]]
    constexpr T length() const noexcept {
        return std::hypot(x, y, z);
    }
    [[using gnu: pure, always_inline, leaf]]
    constexpr T dot(Vec3 obj) const noexcept {
        return x * obj.x + y * obj.y + z * obj.z;
    }
    [[using gnu: pure, always_inline]]
    constexpr Vec3 cross(Vec3 obj) const noexcept {
        return Vec3{y * obj.z - z * obj.y, z * obj.x - x * obj.z, x * obj.y - y * obj.x};
    }
    [[using gnu: pure, always_inline]]
    constexpr T distanceTo(Vec3 obj) const noexcept {
        return std::hypot(x - obj.x, y - obj.y, z - obj.z);
    }

    [[using gnu: pure, always_inline]]
    constexpr Vec3
    operator-() const noexcept {
        return Vec3{-x, -y, -z};
    }

    [[using gnu: pure, always_inline, leaf]]
    constexpr bool
    operator==(Vec3 other) {
        return x == other.x && y == other.y && z == other.z;
    }

    [[using gnu: pure, always_inline]]
    constexpr Vec3
    operator+(Vec3 other) const noexcept {
        return Vec3{x + other.x, y + other.y, z + other.z};
    }
    [[using gnu: always_inline, leaf]]
    constexpr Vec3&
    operator+=(Vec3 other) noexcept {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr Vec3
    operator-(Vec3 other) const noexcept {
        return Vec3{x - other.x, y - other.y, z - other.z};
    }
    [[using gnu: always_inline, leaf]]
    constexpr Vec3&
    operator-=(Vec3 other) noexcept {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    template <Arithmetic U>
    [[using gnu: pure, always_inline]]
    constexpr Vec3
    operator*(U factor) const noexcept {
        return Vec3{
            x * static_cast<T>(factor), y * static_cast<T>(factor), z * static_cast<T>(factor)};
    }
    template <Arithmetic U>
    [[using gnu: always_inline, leaf]]
    constexpr Vec3&
    operator*=(U factor) noexcept {
        x *= static_cast<T>(factor);
        y *= static_cast<T>(factor);
        z *= static_cast<T>(factor);
        return *this;
    }

    template <Arithmetic U>
    [[using gnu: pure, always_inline]]
    constexpr Vec3
    operator/(U den) const noexcept {
        return Vec3{x / static_cast<T>(den), y / static_cast<T>(den), z / static_cast<T>(den)};
    }
    template <Arithmetic U>
    [[using gnu: always_inline, leaf]]
    constexpr Vec3&
    operator/=(U den) noexcept {
        x /= static_cast<T>(den);
        y /= static_cast<T>(den);
        z /= static_cast<T>(den);
        return *this;
    }

    T& x{first_};
    T& y{second_};
    T& z{third_};
};

using Vec3f = Vec3<float>;
using Vec3d = Vec3<double>;

template <std::floating_point T, Arithmetic U>
[[using gnu: const, always_inline]]
constexpr inline Vec3<T>
operator*(U factor, Vec3<T> obj) noexcept {
    return obj * factor;
}

template <std::floating_point T>
[[using gnu: const, always_inline]]
constexpr inline Vec3<T> normalize(Vec3<T> obj) noexcept {
    T length = obj.length();
    return Vec3<T>{obj.x / length, obj.y / length, obj.z / length};
}

} // namespace math
} // namespace boyle

namespace std {

template <floating_point T>
[[using gnu: pure, always_inline]]
constexpr inline T hypot(boyle::math::Vec3<T> obj) noexcept {
    return std::hypot(obj.x, obj.y, obj.z);
}

} // namespace std
