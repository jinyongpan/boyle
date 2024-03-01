/**
 * @file utils.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-07-23
 *
 * @copyright Copyright (c) 2023 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <format>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#include "spdlog/spdlog.h"

#include "math/concepts.hpp"
#include "math/duplet.hpp"
#include "math/triplet.hpp"
#include "math/vec2.hpp"
#include "math/vec3.hpp"

namespace boyle {
namespace math {

constexpr double kEpsilon{1e-8};

template <Arithmetic T>
[[using gnu: const, always_inline, leaf]]
constexpr inline T pow(T x, std::size_t n) noexcept {
    T result{1.0};
    for (std::size_t i{0}; i < n; ++i) {
        result *= x;
    }
    return result;
}

template <std::floating_point T>
[[using gnu: const, always_inline, leaf]]
constexpr inline bool inRange(T value, T start, T end) noexcept {
    if (start > end) {
        start -= end;
        end += start;
        start = end - start;
    }
    return value > start && value < end;
}

template <GeneralArithmetic T>
[[using gnu: const, flatten, leaf]] [[nodiscard]]
inline std::vector<T> linspace(T start, T end, std::size_t num, bool endpoint = true) noexcept {
    if (num == 0) {
        return std::vector<T>(num);
    } else if (num == 1) {
        return std::vector<T>(1, start);
    }
    std::vector<T> result;
    result.reserve(num);
    T step;
    if (endpoint) {
        step = (end - start) / static_cast<double>(num - 1);
    } else {
        step = (end - start) / static_cast<double>(num);
        end -= step;
    }
    for (T x{start}; num; x += step, --num) {
        result.push_back(x);
    }
    return result;
}

template <GeneralArithmetic T, std::floating_point U>
[[using gnu: const, always_inline, leaf, hot]] [[nodiscard]]
constexpr inline T lerp(T start, T end, U ratio) noexcept {
    return (1 - ratio) * start + ratio * end;
}

template <GeneralArithmeticIterator ForwardIt, std::floating_point U = double>
[[using gnu: pure, flatten, leaf, hot]]
inline bool hasDuplicates(ForwardIt first, ForwardIt last, U tol = kEpsilon) noexcept {
    using T = typename std::iterator_traits<ForwardIt>::value_type;
    if (last - first < 2) {
        spdlog::warn(
            "Invalid argument detected! Difference between first and last must be larger than 1: "
            "last - first = {0:d}.",
            last - first
        );
        return false;
    }
    if constexpr (Arithmetic<T>) {
        std::vector<T> sorted{first, last};
        std::sort(sorted.begin(), sorted.end());
        for (typename std::vector<T>::const_iterator it{sorted.cbegin()}; it != sorted.cend() - 1;
             ++it) {
            if (std::abs(*(it + 1) - *it) < tol) {
                return true;
            }
        }
    } else if constexpr (VecArithmetic<T>) {
        for (ForwardIt it{first}; it != last - 1; ++it) {
            if (std::hypot((*(it + 1) - *it)) < tol) {
                return true;
            }
        }
    }
    return false;
}

template <GeneralArithmeticIterator ForwardIt, std::floating_point U = double>
[[using gnu: pure, flatten, leaf, hot]]
inline ForwardIt closetUpperElement(
    ForwardIt first, ForwardIt last, typename std::iterator_traits<ForwardIt>::value_type element,
    U tol = kEpsilon
) noexcept {
    using T = typename std::iterator_traits<ForwardIt>::value_type;

    if constexpr (Arithmetic<T>) {
        if (last - first < 2) {
            spdlog::warn(
                "Invalid argument detected! Difference between first and last must be larger than "
                "1: last - first = {0:d}.",
                last - first
            );
            return element < *first ? first : last;
        }
        if (std::abs(element - *first) < tol) {
            return first + 1;
        } else if (std::abs(element - *(last - 1)) < tol) {
            return last - 1;
        } else {
            return std::upper_bound(first, last, element, [tol](T lhs, T rhs) -> bool {
                return lhs - rhs < -tol;
            });
        }
    } else if constexpr (VecArithmetic<T>) {
        if (last - first < 2) {
            spdlog::error(
                "Invalid argument detected! Difference between first and last must be larger than "
                "1: last - first = {0:d}.",
                last - first
            );
            return ForwardIt{};
        }
        const std::size_t size = last - first;
        std::size_t pos;
        U min_distance = std::numeric_limits<U>::max();
        for (ForwardIt it{first}; it != last; ++it) {
            const U distance = element.distanceTo(*it);
            if (distance < min_distance) {
                pos = it - first;
                min_distance = distance;
            }
        }
        T diff, r;
        U inner_prod;
        if (pos == 0) {
            diff = *(first + 1) - *first;
            r = element - *first;
            inner_prod = diff.dot(r);
            if (inner_prod < -tol) {
                return first;
            } else {
                return first + 1;
            }
        } else if (pos == size - 1) {
            diff = *(last - 2) - *(last - 1);
            r = element - *(last - 1);
            inner_prod = diff.dot(r);
            if (inner_prod < -tol) {
                return last;
            } else {
                return last - 1;
            }
        } else {
            ForwardIt it = first + pos;
            diff = *(it + 1) - *it;
            r = element - *it;
            inner_prod = diff.dot(r);
            if (inner_prod < -tol) {
                return it;
            } else {
                return it + 1;
            }
        }
    }
}

template <typename T0, typename T1>
[[using gnu: pure, flatten, leaf]] [[nodiscard]]
inline std::vector<Duplet<T0, T1>> stack(
    const std::vector<T0>& xs, const std::vector<T1>& ys
) noexcept {
    if (xs.size() != ys.size()) {
        spdlog::error(
            "Invalid arguments detected! xs, ys must share the same size: xs.size() = {0:d} while "
            "ys.size() = {1:d}",
            xs.size(), ys.size()
        );
        return std::vector<Duplet<T0, T1>>{};
    }
    const std::size_t size{xs.size()};
    std::vector<Duplet<T0, T1>> result;
    result.reserve(size);
    for (std::size_t i{0}; i < size; ++i) {
        result.emplace_back(xs[i], ys[i]);
    }
    return result;
}

template <typename T0, typename T1, typename T2>
[[using gnu: pure, flatten, leaf]] [[nodiscard]]
inline std::vector<Triplet<T0, T1, T2>> stack(
    const std::vector<T0>& xs, const std::vector<T1>& ys, const std::vector<T2>& zs
) noexcept {
    if (xs.size() != ys.size() || ys.size() != zs.size()) {
        spdlog::error(
            "Invalid arguments detected! xs, ys must share the same size: xs.size() = {0:d} while "
            "ys.size() = {1:d}",
            xs.size(), ys.size()
        );
        return std::vector<Triplet<T0, T1, T2>>{};
    }
    const std::size_t size{xs.size()};
    std::vector<Triplet<T0, T1, T2>> result;
    result.reserve(size);
    for (std::size_t i{0}; i < size; ++i) {
        result.emplace_back(xs[i], ys[i], zs[i]);
    }
    return result;
}

} // namespace math
} // namespace boyle
