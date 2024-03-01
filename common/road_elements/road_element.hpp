/**
 * @file road_element.h
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-20
 *
 * @copyright Copyright (c) 2023 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include <cstdint>

#include "common/utils/macros.hpp"

namespace boyle {
namespace common {

class RoadElement {
  public:
    explicit RoadElement(std::uint64_t id) noexcept : id_(id) {}
    ENABLE_IMPLICIT_CONSTRUCTORS(RoadElement);
    virtual ~RoadElement() noexcept = 0;
    const std::uint64_t& id() const noexcept { return id_; }

  protected:
    std::uint64_t id_;
};

inline RoadElement::~RoadElement() noexcept = default;

} // namespace common
} // namespace boyle
