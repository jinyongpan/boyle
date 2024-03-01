/**
 * @file triplet.hpp
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
#include <tuple>
#include <utility>

#include "boost/serialization/access.hpp"

#define DECLARE_TRIPLET(ClassName, VAR_0, VAR_1, VAR_2)                                         \
    template <typename T0, typename T1 = T0, typename T2 = T0>                                  \
    class ClassName final : public ::boyle::math::Triplet<T0, T1, T2> {                       \
        using ::boyle::math::Triplet<T0, T1, T2>::first_;                                     \
        using ::boyle::math::Triplet<T0, T1, T2>::second_;                                    \
        using ::boyle::math::Triplet<T0, T1, T2>::third_;                                     \
                                                                                                \
      public:                                                                                   \
        [[using gnu: always_inline]] ClassName() noexcept                                       \
            : ::boyle::math::Triplet<T0, T1, T2>{} {}                                         \
        [[using gnu: always_inline]] constexpr ClassName(T0 cv) noexcept                        \
            requires std::same_as<T0, T1> && std::same_as<T0, T2>                               \
            : ::boyle::math::Triplet<T0>{cv} {}                                               \
        [[using gnu: always_inline]] constexpr ClassName(                                       \
            T0 c##VAR_0, T1 c##VAR_1, T2 c##VAR_2                                               \
        ) noexcept                                                                              \
            : ::boyle::math::Triplet<T0, T1, T2>{                                             \
                  std::move(c##VAR_0), std::move(c##VAR_1), std::move(c##VAR_2)} {}             \
        [[using gnu: always_inline]] constexpr ClassName(const ClassName& other) noexcept       \
            : ::boyle::math::Triplet<T0, T1, T2>{other} {}                                    \
        [[using gnu: always_inline]]                                                            \
        constexpr ClassName&                                                                    \
        operator=(const ClassName& other) noexcept {                                            \
            ::boyle::math::Triplet<T0, T1, T2>::operator=(other);                             \
            return *this;                                                                       \
        }                                                                                       \
        ~ClassName() noexcept override {}                                                       \
                                                                                                \
        [[using gnu: always_inline]] constexpr ClassName(std::tuple<T0, T1, T2> tuple) noexcept \
            : ::boyle::math::Triplet<T0, T1, T2>{std::move(tuple)} {}                         \
                                                                                                \
        template <                                                                              \
            std::convertible_to<T0> U0, std::convertible_to<T1> U1, std::convertible_to<T2> U2> \
        [[using gnu: pure, always_inline]] constexpr operator ClassName<U0, U1, U2>(            \
        ) const noexcept {                                                                      \
            return ClassName<T0, T1, T2>{                                                       \
                static_cast<U0>(first_), static_cast<U1>(second_), static_cast<U2>(third_)};    \
        }                                                                                       \
                                                                                                \
        T0& VAR_0{first_};                                                                      \
        T1& VAR_1{second_};                                                                     \
        T2& VAR_2{third_};                                                                      \
    }

namespace boyle {
namespace math {

template <typename T0, typename T1 = T0, typename T2 = T0>
class Triplet {
    template <typename U0, typename U1, typename U2>
    friend inline std::ostream& std::operator<<(
        std::ostream& os, const Triplet<U0, U1, U2>& obj
    ) noexcept;
    friend class boost::serialization::access;

  public:
    [[using gnu: always_inline]] Triplet() noexcept = default;
    [[using gnu: always_inline]] constexpr Triplet(T0 cv) noexcept
        requires std::same_as<T0, T1> && std::same_as<T0, T2>
        : first_{cv}, second_{cv}, third_{cv} {}
    [[using gnu: always_inline]] constexpr Triplet(T0 cv0, T1 cv1, T2 cv2) noexcept
        : first_{std::move(cv0)}, second_{std::move(cv1)}, third_{std::move(cv2)} {}
    [[using gnu: always_inline]] constexpr Triplet(const Triplet& other) noexcept = default;
    [[using gnu: always_inline]]
    constexpr Triplet&
    operator=(const Triplet& other) noexcept = default;
    virtual ~Triplet() noexcept = 0;

    [[using gnu: always_inline]] constexpr Triplet(std::tuple<T0, T1, T2> tuple) noexcept
        : first_{std::move(std::get<0>(tuple))}, second_{std::move(std::get<1>(tuple))},
          third_{std::move(std::get<2>(tuple))} {}

    [[using gnu: pure, always_inline]] constexpr operator std::tuple<T0, T1, T2>() const noexcept {
        return std::make_tuple<T0, T1, T2>(first_, second_, third_);
    }

  protected:
    template <typename Archive>
    [[using gnu: always_inline]]
    void serialize(Archive& ar, const unsigned int version) noexcept {
        ar& first_;
        ar& second_;
        ar& third_;
        return;
    }

    T0 first_;
    T1 second_;
    T2 third_;
};

template <typename T0, typename T1, typename T2>
inline Triplet<T0, T1, T2>::~Triplet() noexcept = default;

DECLARE_TRIPLET(SlzTriplet, s, l, z);

using SlzTripletf = SlzTriplet<float>;
using SlzTripletd = SlzTriplet<double>;

} // namespace math
} // namespace boyle

namespace std {

template <typename T0, typename T1, typename T2>
[[using gnu: always_inline]]
inline ostream&
operator<<(ostream& os, const boyle::math::Triplet<T0, T1, T2>& obj) noexcept {
    os << "(" << obj.first_ << ", " << obj.second_ << ", " << obj.third_ << ")";
    return os;
}

} // namespace std
