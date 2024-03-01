/**
 * @file vec2.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-07-17
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
#include "math/duplet.hpp"

namespace boyle {
namespace math {

template <std::floating_point T>
class Vec2 final : public Duplet<T> {
    using Duplet<T>::first_;
    using Duplet<T>::second_;

  public:
    using value_type = T;

    [[using gnu: always_inline]] Vec2() noexcept : Duplet<T>{} {}
    [[using gnu: always_inline]] constexpr Vec2(T cv) noexcept : Duplet<T>{cv, cv} {}
    [[using gnu: always_inline]] constexpr Vec2(T cx, T cy) noexcept : Duplet<T>{cx, cy} {}
    [[using gnu: always_inline]] constexpr Vec2(const Vec2& other) noexcept : Duplet<T>{other} {}
    [[using gnu: always_inline]]
    constexpr Vec2&
    operator=(const Vec2& other) noexcept {
        Duplet<T>::operator=(other);
        return *this;
    };
    ~Vec2() noexcept override {}

    template <std::floating_point U>
    [[using gnu: always_inline]] constexpr Vec2(const std::pair<U, U>& other) noexcept
        : Duplet<T>{other} {}

    template <std::floating_point U>
    [[using gnu: pure, always_inline]] constexpr operator Vec2<U>() const noexcept {
        return Vec2<U>{static_cast<U>(x), static_cast<U>(y)};
    }

    [[using gnu: pure, always_inline]]
    constexpr T length() const noexcept {
        return std::hypot(x, y);
    }
    [[using gnu: pure, always_inline]]
    constexpr T lengthSqr() const noexcept {
        return x * x + y * y;
    }
    [[using gnu: pure, always_inline]]
    constexpr T angle() const noexcept {
        return std::atan2(y, x);
    }
    [[using gnu: pure, always_inline, leaf]]
    constexpr T dot(Vec2 obj) const noexcept {
        return x * obj.x + y * obj.y;
    }
    [[using gnu: pure, always_inline, leaf]]
    constexpr T cross(Vec2 obj) const noexcept {
        return x * obj.y - y * obj.x;
    }
    [[using gnu: pure, always_inline]]
    constexpr T distanceTo(Vec2 obj) const noexcept {
        return std::hypot(x - obj.x, y - obj.y);
    }

    template <std::floating_point U>
    [[using gnu: pure, always_inline]]
    constexpr Vec2 rotate(U radian) const noexcept {
        return Vec2{
            x * std::cos(radian) - y * std::sin(radian),
            x * std::sin(radian) + y * std::cos(radian)};
    }
    template <std::floating_point U>
    [[using gnu: pure, always_inline]]
    constexpr Vec2& selfRotate(U radian) noexcept {
        *this = Vec2{
            x * std::cos(radian) - y * std::sin(radian),
            x * std::sin(radian) + y * std::cos(radian)};
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr Vec2 rotateHalfPi() const noexcept {
        return Vec2{-y, x};
    }
    [[using gnu: always_inline]]
    constexpr Vec2& selfRotateHalfPi() noexcept {
        *this = Vec2{-y, x};
        return *this;
    }

    [[using gnu: pure, always_inline, leaf]]
    constexpr bool
    operator==(Vec2 other) {
        return x == other.x && y == other.y;
    }

    [[using gnu: pure, always_inline]]
    constexpr Vec2
    operator-() const noexcept {
        return Vec2{-x, -y};
    }

    [[using gnu: pure, always_inline]]
    constexpr Vec2
    operator+(Vec2 other) const noexcept {
        return Vec2{x + other.x, y + other.y};
    }
    [[using gnu: always_inline, leaf]]
    constexpr Vec2&
    operator+=(Vec2 other) noexcept {
        x += other.x;
        y += other.y;
        return *this;
    }

    [[using gnu: pure, always_inline]]
    constexpr Vec2
    operator-(Vec2 other) const noexcept {
        return Vec2{x - other.x, y - other.y};
    }
    [[using gnu: always_inline, leaf]]
    constexpr Vec2&
    operator-=(Vec2 other) noexcept {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    template <Arithmetic U>
    [[using gnu: pure, always_inline]]
    constexpr Vec2
    operator*(U factor) const noexcept {
        return Vec2{x * static_cast<T>(factor), y * static_cast<T>(factor)};
    }
    template <Arithmetic U>
    [[using gnu: always_inline, leaf]]
    constexpr Vec2&
    operator*=(U factor) noexcept {
        x *= static_cast<T>(factor);
        y *= static_cast<T>(factor);
        return *this;
    }

    template <Arithmetic U>
    [[using gnu: pure, always_inline]]
    constexpr Vec2
    operator/(U den) const noexcept {
        return Vec2{x / static_cast<T>(den), y / static_cast<T>(den)};
    }
    template <Arithmetic U>
    [[using gnu: always_inline, leaf]]
    constexpr Vec2&
    operator/=(U den) noexcept {
        x /= static_cast<T>(den);
        y /= static_cast<T>(den);
        return *this;
    }

    T& x{first_};
    T& y{second_};
};

using Vec2f = Vec2<float>;
using Vec2d = Vec2<double>;

template <std::floating_point T, Arithmetic U>
[[using gnu: const, always_inline]]
constexpr inline Vec2<T>
operator*(U factor, Vec2<T> obj) noexcept {
    return obj * factor;
}

template <std::floating_point T>
[[using gnu: const, always_inline]]
constexpr inline Vec2<T> normalize(Vec2<T> obj) noexcept {
    const T length = obj.length();
    return Vec2<T>{obj.x / length, obj.y / length};
}

} // namespace math
} // namespace boyle

namespace std {

template <floating_point T>
[[using gnu: const, always_inline]]
constexpr inline T hypot(boyle::math::Vec2<T> obj) noexcept {
    return std::hypot(obj.x, obj.y);
}

template <floating_point T>
[[using gnu: const, always_inline]]
constexpr inline T atan2(boyle::math::Vec2<T> obj) noexcept {
    return std::atan2(obj.y, obj.x);
}

} // namespace std
