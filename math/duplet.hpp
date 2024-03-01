/**
 * @file duplet.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-12-16
 *
 * @copyright Copyright (c) 2023 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include <concepts>
#include <ostream>
#include <utility>

#include "boost/serialization/access.hpp"

#define DECLARE_DUPLET(ClassName, VAR_0, VAR_1)                                                    \
    template <typename T0, typename T1 = T0>                                                       \
    class ClassName final : public ::boyle::math::Duplet<T0, T1> {                               \
        using ::boyle::math::Duplet<T0, T1>::first_;                                             \
        using ::boyle::math::Duplet<T0, T1>::second_;                                            \
                                                                                                   \
      public:                                                                                      \
        [[using gnu: always_inline]] ClassName() noexcept : ::boyle::math::Duplet<T0, T1>{} {}   \
        [[using gnu: always_inline]] constexpr ClassName(T0 cv) noexcept                           \
            requires std::same_as<T0, T1>                                                          \
            : ::boyle::math::Duplet<T0>{cv} {}                                                   \
        [[using gnu: always_inline]] constexpr ClassName(T0 c##VAR_0, T1 c##VAR_1) noexcept        \
            : ::boyle::math::Duplet<T0, T1>{std::move(c##VAR_0), std::move(c##VAR_1)} {}         \
        [[using gnu: always_inline]] constexpr ClassName(const ClassName& other) noexcept          \
            : ::boyle::math::Duplet<T0, T1>{other} {}                                            \
        [[using gnu: always_inline]]                                                               \
        constexpr ClassName&                                                                       \
        operator=(const ClassName& other) noexcept {                                               \
            ::boyle::math::Duplet<T0, T1>::operator=(other);                                     \
            return *this;                                                                          \
        }                                                                                          \
        ~ClassName() noexcept override {}                                                          \
                                                                                                   \
        [[using gnu: always_inline]] constexpr ClassName(std::pair<T0, T1> pair) noexcept          \
            : ::boyle::math::Duplet<T0, T1>{std::move(pair)} {}                                  \
                                                                                                   \
        template <std::convertible_to<T0> U0, std::convertible_to<T1> U1>                          \
        [[using gnu: pure, always_inline]] constexpr operator ClassName<U0, U1>() const noexcept { \
            return ClassName<U0, U1>{static_cast<U0>(first_), static_cast<U1>(second_)};           \
        }                                                                                          \
                                                                                                   \
        T0& VAR_0{first_};                                                                         \
        T1& VAR_1{second_};                                                                        \
    }

namespace boyle {
namespace math {

template <typename T0, typename T1 = T0>
class Duplet {
    template <typename U0, typename U1>
    friend inline std::ostream& std::operator<<(
        std::ostream& os, const Duplet<U0, U1>& obj
    ) noexcept;
    friend class boost::serialization::access;

  public:
    [[using gnu: always_inline]] Duplet() noexcept = default;
    [[using gnu: always_inline]] constexpr Duplet(T0 cv) noexcept
        requires std::same_as<T0, T1>
        : first_{cv}, second_{cv} {}
    [[using gnu: always_inline]] constexpr Duplet(T0 cv0, T1 cv1) noexcept
        : first_{std::move(cv0)}, second_{std::move(cv1)} {}
    [[using gnu: always_inline]] constexpr Duplet(const Duplet& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr Duplet&
    operator=(const Duplet& other) noexcept = default;
    virtual ~Duplet() noexcept = 0;

    [[using gnu: always_inline]] constexpr Duplet(std::pair<T0, T1> pair) noexcept
        : first_{std::move(pair.first)}, second_{std::move(pair.second)} {}

    [[using gnu: pure, always_inline]] constexpr operator std::pair<T0, T1>() const noexcept {
        return std::make_pair<T0, T1>(first_, second_);
    }

  protected:
    template <typename Archive>
    [[using gnu: always_inline]]
    void serialize(Archive& ar, const unsigned int version) noexcept {
        ar& first_;
        ar& second_;
        return;
    }

    T0 first_;
    T1 second_;
};

template <typename T0, typename T1>
inline Duplet<T0, T1>::~Duplet() noexcept = default;

DECLARE_DUPLET(SlDuplet, s, l);

DECLARE_DUPLET(StDuplet, s, t);

DECLARE_DUPLET(SxDuplet, s, x);

DECLARE_DUPLET(SyDuplet, s, y);

using SlDupletf = SlDuplet<float>;
using SlDupletd = SlDuplet<double>;

using StDupletf = StDuplet<float>;
using StDupletd = StDuplet<double>;

using SxDupletf = SxDuplet<float>;
using SxDupletd = SxDuplet<double>;

using SyDupletf = SyDuplet<float>;
using SyDupletd = SyDuplet<double>;

} // namespace math
} // namespace boyle

namespace std {

template <typename T0, typename T1>
[[using gnu: always_inline]]
inline ostream&
operator<<(ostream& os, const boyle::math::Duplet<T0, T1>& obj) noexcept {
    os << "(" << obj.first_ << ", " << obj.second_ << ")";
    return os;
}

} // namespace std
