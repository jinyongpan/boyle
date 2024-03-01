/**
 * @file concepts.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-11-15
 *
 * @copyright Copyright (c) 2023 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>
#include <iterator>
#include <type_traits>

namespace boyle {
namespace math {

template <typename T>
concept Arithmetic = std::is_arithmetic_v<T>;

template <typename T>
concept VecArithmetic = requires(T a, T b) {
    { a + b } -> std::same_as<T>;
    { a - b } -> std::same_as<T>;
    { a.dot(b) } -> std::floating_point;
    { a * 1.0 } -> std::same_as<T>;
    { 1.0 * a } -> std::same_as<T>;
    { a / 1.0 } -> std::same_as<T>;
    { a += b } -> std::same_as<T&>;
    { a -= b } -> std::same_as<T&>;
    { a.distanceTo(b) } -> std::floating_point;
};

template <typename T>
concept GeneralArithmetic = Arithmetic<T> || VecArithmetic<T>;

template <typename T>
concept ArithmeticIterator = Arithmetic<typename std::iterator_traits<T>::value_type>;

template <typename T>
concept VecArithmeticIterator = VecArithmetic<typename std::iterator_traits<T>::value_type>;

template <typename T>
concept GeneralArithmeticIterator = ArithmeticIterator<T> || VecArithmeticIterator<T>;

template <typename T>
concept FloatingPointIterator = std::floating_point<typename std::iterator_traits<T>::value_type>;

template <typename T, typename U>
concept SameValueTypeIterator = std::same_as<typename std::iterator_traits<T>::value_type, U>;

} // namespace math
} // namespace boyle
