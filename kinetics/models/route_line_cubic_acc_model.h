/**
 * @file route_line_cubic_acc_model.h
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-15
 *
 * @copyright Copyright (c) 2023 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include <numeric>
#include <vector>

#include "common/utils/macros.hpp"
#include "kinetics/fence1.hpp"
#include "kinetics/motion1.hpp"
#include "math/qp_solvers/osqp_solver.hpp"
#include "math/qp_solvers/qp_problem.hpp"

namespace boyle {
namespace kinetics {

class [[nodiscard]] RouteLineCubicAccModel final {
  public:
    struct [[nodiscard]] Result final {
        Motion1d motion;
        ::boyle::math::OsqpInfo info;
    };

    explicit RouteLineCubicAccModel(std::vector<double> sample_ts);
    explicit RouteLineCubicAccModel(
        std::vector<double> sample_ts, ::boyle::math::OsqpSettings settings
    );
    DISABLE_COPY_AND_MOVE(RouteLineCubicAccModel);
    ~RouteLineCubicAccModel() noexcept = default;
    const std::size_t& num_samples() const noexcept;
    const ::boyle::math::QpProblem<double, int>& qp_problem() const noexcept;
    const ::boyle::math::OsqpSettings& settings() const noexcept;
    void setSoftFences(const std::vector<SoftFence1d>& soft_fences) noexcept;
    void setHardFences(const std::vector<HardFence1d>& hard_fences) noexcept;
    void setVelocityRange(double lower_bound, double upper_bound) noexcept;
    void setAccelRange(double lower_bound, double upper_bound) noexcept;
    void setInitialState(
        double s0 = std::numeric_limits<double>::quiet_NaN(),
        double v0 = std::numeric_limits<double>::quiet_NaN(),
        double a0 = std::numeric_limits<double>::quiet_NaN()
    ) noexcept;
    void setFinalState(
        double sf = std::numeric_limits<double>::quiet_NaN(),
        double vf = std::numeric_limits<double>::quiet_NaN(),
        double af = std::numeric_limits<double>::quiet_NaN()
    ) noexcept;
    void setVelocityCost(double target_velocity, double velocity_weight) noexcept;
    void setAccelCost(double accel_weight) noexcept;
    void setJerkCost(double jerk_weight) noexcept;
    Result solve() const noexcept;
    void clear() noexcept;

  private:
    static std::array<double, 4> prorationCoeffs(double ratio, double scale = 1.0) noexcept;
    void setIntegrationRelation() noexcept;
    std::size_t sIndex(std::size_t t_index) const noexcept;
    std::size_t vIndex(std::size_t t_index) const noexcept;
    std::size_t num_samples_;
    double time_scale_;
    ::boyle::math::QpProblem<double, int> qp_problem_;
    ::boyle::math::OsqpSettings settings_;
    std::vector<double> sample_ts_;
    std::vector<double> hs_;
    std::vector<double> h2s_;
    std::vector<double> reciprocal_hs_;
    std::vector<double> reciprocal_h2s_;
};

} // namespace kinetics
} // namespace boyle
