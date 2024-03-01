/**
 * @file route_line_cubic_acc_model.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-11-23
 *
 * @copyright Copyright (c) 2023 Boyle Development Team.
 *            All rights reserved.
 *
 */

#include "route_line_cubic_acc_model.h"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

#include "spdlog/spdlog.h"

namespace boyle {
namespace kinetics {

std::array<double, 4> RouteLineCubicAccModel::prorationCoeffs(double ratio, double scale) noexcept {
    const double ratio2 = ratio * ratio;
    const double ratio3 = ratio2 * ratio;
    return std::array<double, 4>{
        1.0 - ratio2 * 3.0 + ratio3 * 2.0, ratio2 * 3.0 - ratio3 * 2.0,
        (ratio - ratio2 * 2.0 + ratio3) * scale, -(ratio2 - ratio3) * scale};
}

RouteLineCubicAccModel::RouteLineCubicAccModel(std::vector<double> sample_ts) {
    if (sample_ts.size() < 2) {
        std::string error_msg = std::format(
            "Invalid argument error detected: size of sample_ts must larger than 2: "
            "sample_ts.size() = {0:d}.",
            sample_ts.size()
        );
        throw std::invalid_argument(std::move(error_msg));
    }
    num_samples_ = sample_ts.size();
    qp_problem_.resize(num_samples_ * 2, num_samples_ * 5);
    osqp_set_default_settings(&settings_);
    settings_.scaling = 0;
    sample_ts_ = std::move(sample_ts);
    time_scale_ = sample_ts_.back() - sample_ts_.front();
    hs_.reserve(num_samples_ - 1);
    h2s_.reserve(num_samples_ - 1);
    reciprocal_hs_.reserve(num_samples_ - 1);
    reciprocal_h2s_.reserve(num_samples_ - 1);
    for (std::vector<double>::const_iterator it = sample_ts_.cbegin() + 1; it != sample_ts_.cend();
         ++it) {
        const double diff = *it - *(it - 1);
        const double reciprocal_diff = 1.0 / diff;
        hs_.push_back(diff);
        h2s_.push_back(diff * diff);
        reciprocal_hs_.push_back(reciprocal_diff);
        reciprocal_h2s_.push_back(reciprocal_diff * reciprocal_diff);
    }
    setIntegrationRelation();
    setInitialState();
    setFinalState();
    for (std::size_t i = 1; i < num_samples_ - 1; ++i) {
        qp_problem_.updateConstrainTerm(
            sIndex(i), {{sIndex(i), 1.0}}, std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::max()
        );
        qp_problem_.updateConstrainTerm(
            vIndex(i), {{vIndex(i), 1.0}}, std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::max()
        );
    }
}

RouteLineCubicAccModel::RouteLineCubicAccModel(
    std::vector<double> sample_ts, ::boyle::math::OsqpSettings settings
)
    : RouteLineCubicAccModel(std::move(sample_ts)) {
    settings_ = settings;
}

const std::size_t& RouteLineCubicAccModel::num_samples() const noexcept { return num_samples_; }

const ::boyle::math::QpProblem<double, int>& RouteLineCubicAccModel::qp_problem() const noexcept {
    return qp_problem_;
}

const ::boyle::math::OsqpSettings& RouteLineCubicAccModel::settings() const noexcept {
    return settings_;
}

void RouteLineCubicAccModel::setIntegrationRelation() noexcept {
    for (std::size_t i{1}; i < num_samples_ - 1; ++i) {
        const std::unordered_map<int, double> constrain_vec{
            {sIndex(i - 1), reciprocal_h2s_[i - 1] * 3.0},
            {sIndex(i), -(reciprocal_h2s_[i - 1] - reciprocal_h2s_[i]) * 3.0},
            {sIndex(i + 1), -(reciprocal_h2s_[i] * 3.0)},
            {vIndex(i - 1), reciprocal_hs_[i - 1]},
            {vIndex(i), (reciprocal_hs_[i - 1] + reciprocal_hs_[i]) * 2.0},
            {vIndex(i + 1), reciprocal_hs_[i]}};
        qp_problem_.updateConstrainTerm(
            num_samples_ * 3 + i, constrain_vec, -boyle::math::kEpsilon, ::boyle::math::kEpsilon
        );
    }
    return;
}

void RouteLineCubicAccModel::setSoftFences(const std::vector<SoftFence1d>& soft_fences) noexcept {
    for (const SoftFence1d& soft_fence : soft_fences) {
        const std::size_t istart =
            ::boyle::math::closetUpperElement(
                sample_ts_.cbegin(), sample_ts_.cend(), soft_fence.bound_line.minT()
            ) -
            sample_ts_.cbegin();
        const std::size_t iend =
            ::boyle::math::closetUpperElement(
                sample_ts_.cbegin(), sample_ts_.cend(), soft_fence.bound_line.maxT()
            ) -
            sample_ts_.cbegin();
        if (istart == num_samples_) {
            spdlog::warn(
                "Invalid argument issue detected! The minT() of hard_fence.bound_line should be "
                "less than sample_ts_.back(): hard_fence.bound_line.minT() = {0:f} while "
                "sample_ts_.back() = {1:f}.",
                soft_fence.bound_line.minT(), sample_ts_.back()
            );
            break;
        }
        if (iend == 0) {
            spdlog::warn(
                "Invalid argument issue detected! The maxT() of hard_fence.bound_line should be "
                "larger than sample_ts_.front(): hard_fence.bound_line.maxT() = {0:f} while "
                "sample_ts_.front() = {1:f}.",
                soft_fence.bound_line.maxT(), sample_ts_.front()
            );
            break;
        }

        double factor;

        if (soft_fence.actio == ::boyle::common::Actio::PUSHING) {
            if (istart == iend) {
                const double h{sample_ts_[istart] - sample_ts_[istart - 1]};
                double ratio{(soft_fence.bound_line.minT() - sample_ts_[istart - 1]) / h};
                std::array<double, 4> proration_coeffs = prorationCoeffs(ratio, h);
                std::unordered_map<int, double> constrain_vec{
                    {sIndex(istart - 1), -proration_coeffs[0]},
                    {sIndex(istart), -proration_coeffs[1]},
                    {vIndex(istart - 1), -proration_coeffs[2]},
                    {vIndex(istart), -proration_coeffs[3]}};
                factor = (soft_fence.bound_line.maxT() - soft_fence.bound_line.minT()) /
                         time_scale_ * 0.5;
                qp_problem_.addClampCostTerm(
                    constrain_vec, -soft_fence.bound_line.ys().front(),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
                ratio = (soft_fence.bound_line.maxT() - sample_ts_[iend - 1]) / h;
                proration_coeffs = prorationCoeffs(ratio, h);
                constrain_vec = {
                    {sIndex(iend - 1), -proration_coeffs[0]},
                    {sIndex(iend), -proration_coeffs[1]},
                    {vIndex(iend - 1), -proration_coeffs[2]},
                    {vIndex(iend), -proration_coeffs[3]}};
                qp_problem_.addClampCostTerm(
                    constrain_vec, -soft_fence.bound_line.ys().back(),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
                break;
            }
            if (istart != 0) {
                const double h{sample_ts_[istart] - sample_ts_[istart - 1]};
                const double ratio{(soft_fence.bound_line.minT() - sample_ts_[istart - 1]) / h};
                const std::array<double, 4> proration_coeffs = prorationCoeffs(ratio, h);
                const std::unordered_map<int, double> constrain_vec{
                    {sIndex(istart - 1), -proration_coeffs[0]},
                    {sIndex(istart), -proration_coeffs[1]},
                    {vIndex(istart - 1), -proration_coeffs[2]},
                    {vIndex(istart), -proration_coeffs[3]}};
                factor = (sample_ts_[istart] - soft_fence.bound_line.minT()) / time_scale_;
                qp_problem_.addClampCostTerm(
                    constrain_vec, -soft_fence.bound_line.ys().front(),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
            }
            for (std::size_t i{istart}; i < iend; ++i) {
                if (i == istart) {
                    factor = hs_[i] / time_scale_ * 0.5;
                } else if (i == iend - 1) {
                    factor = hs_[i - 1] / time_scale_ * 0.5;
                } else {
                    factor = (hs_[i - 1] + hs_[i]) / time_scale_ * 0.5;
                }
                qp_problem_.addClampCostTerm(
                    {{sIndex(i), -1.0}}, -soft_fence.bound_line(sample_ts_[i]),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
            }
            if (iend != num_samples_) {
                const double h{sample_ts_[iend] - sample_ts_[iend - 1]};
                const double ratio{(soft_fence.bound_line.maxT() - sample_ts_[iend - 1]) / h};
                const std::array<double, 4> proration_coeffs = prorationCoeffs(ratio, h);
                const std::unordered_map<int, double> constrain_vec{
                    {sIndex(iend - 1), -proration_coeffs[0]},
                    {sIndex(iend), -proration_coeffs[1]},
                    {vIndex(iend - 1), -proration_coeffs[2]},
                    {vIndex(iend), -proration_coeffs[3]}};
                factor = (soft_fence.bound_line.maxT() - sample_ts_[iend - 1]) / time_scale_;
                qp_problem_.addClampCostTerm(
                    constrain_vec, -soft_fence.bound_line.ys().back(),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
            }
        } else {
            if (istart == iend) {
                const double h{sample_ts_[istart] - sample_ts_[istart - 1]};
                double ratio{(soft_fence.bound_line.minT() - sample_ts_[istart - 1]) / h};
                std::array<double, 4> proration_coeffs = prorationCoeffs(ratio, h);
                std::unordered_map<int, double> constrain_vec{
                    {sIndex(istart - 1), proration_coeffs[0]},
                    {sIndex(istart), proration_coeffs[1]},
                    {vIndex(istart - 1), proration_coeffs[2]},
                    {vIndex(istart), proration_coeffs[3]}};
                factor = (soft_fence.bound_line.maxT() - soft_fence.bound_line.minT()) /
                         time_scale_ * 0.5;
                qp_problem_.addClampCostTerm(
                    constrain_vec, soft_fence.bound_line.ys().front(),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
                ratio = (soft_fence.bound_line.maxT() - sample_ts_[iend - 1]) / h;
                proration_coeffs = prorationCoeffs(ratio, h);
                constrain_vec = {
                    {sIndex(iend - 1), proration_coeffs[0]},
                    {sIndex(iend), proration_coeffs[1]},
                    {vIndex(iend - 1), proration_coeffs[2]},
                    {vIndex(iend), proration_coeffs[3]}};
                qp_problem_.addClampCostTerm(
                    constrain_vec, soft_fence.bound_line.ys().back(),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
                break;
            }
            if (istart != 0) {
                const double h{sample_ts_[istart] - sample_ts_[istart - 1]};
                const double ratio{(soft_fence.bound_line.minT() - sample_ts_[istart - 1]) / h};
                const std::array<double, 4> proration_coeffs = prorationCoeffs(ratio, h);
                const std::unordered_map<int, double> constrain_vec{
                    {sIndex(istart - 1), proration_coeffs[0]},
                    {sIndex(istart), proration_coeffs[1]},
                    {vIndex(istart - 1), proration_coeffs[2]},
                    {vIndex(istart), proration_coeffs[3]}};
                factor = (sample_ts_[istart] - soft_fence.bound_line.minT()) / time_scale_;
                qp_problem_.addClampCostTerm(
                    constrain_vec, soft_fence.bound_line.ys().front(),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
            }
            for (std::size_t i{istart}; i < iend; ++i) {
                if (i == istart) {
                    factor = hs_[i] / time_scale_ * 0.5;
                } else if (i == iend - 1) {
                    factor = hs_[i - 1] / time_scale_ * 0.5;
                } else {
                    factor = (hs_[i - 1] + hs_[i]) / time_scale_ * 0.5;
                }
                qp_problem_.addClampCostTerm(
                    {{sIndex(i), 1.0}}, soft_fence.bound_line(sample_ts_[i]),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
            }
            if (iend != num_samples_) {
                const double h{sample_ts_[iend] - sample_ts_[iend - 1]};
                const double ratio{(soft_fence.bound_line.maxT() - sample_ts_[iend - 1]) / h};
                const std::array<double, 4> proration_coeffs = prorationCoeffs(ratio, h);
                const std::unordered_map<int, double> constrain_vec{
                    {sIndex(iend - 1), proration_coeffs[0]},
                    {sIndex(iend), proration_coeffs[1]},
                    {vIndex(iend - 1), proration_coeffs[2]},
                    {vIndex(iend), proration_coeffs[3]}};
                factor = (soft_fence.bound_line.maxT() - sample_ts_[iend - 1]) / time_scale_;
                qp_problem_.addClampCostTerm(
                    constrain_vec, soft_fence.bound_line.ys().back(),
                    soft_fence.linear_weight * factor, soft_fence.quadratic_weight * factor
                );
            }
        }
    }
    return;
}

void RouteLineCubicAccModel::setHardFences(const std::vector<HardFence1d>& hard_fences) noexcept {
    std::vector<double> lower_bound(num_samples_, std::numeric_limits<double>::lowest());
    std::vector<double> upper_bound(num_samples_, std::numeric_limits<double>::max());
    for (const HardFence1d& hard_fence : hard_fences) {
        const std::size_t istart =
            ::boyle::math::closetUpperElement(
                sample_ts_.cbegin(), sample_ts_.cend(), hard_fence.bound_line.minT()
            ) -
            sample_ts_.cbegin();
        const std::size_t iend =
            ::boyle::math::closetUpperElement(
                sample_ts_.cbegin(), sample_ts_.cend(), hard_fence.bound_line.maxT()
            ) -
            sample_ts_.cbegin();
        if (istart == num_samples_) {
            spdlog::warn(
                "Invalid argument issue detected! The minT() of hard_fence.bound_line should be "
                "less than sample_ts_.back(): hard_fence.bound_line.minT() = {0:f} while "
                "sample_ts_.back() = {1:f}.",
                hard_fence.bound_line.minT(), sample_ts_.back()
            );
            break;
        }
        if (iend == 0) {
            spdlog::warn(
                "Invalid argument issue detected! The maxT() of hard_fence.bound_line should be "
                "larger than sample_ts_.front(): hard_fence.bound_line.maxT() = {0:f} while "
                "sample_ts_.front() = {1:f}.",
                hard_fence.bound_line.maxT(), sample_ts_.front()
            );
            break;
        }

        if (hard_fence.actio == ::boyle::common::Actio::PUSHING) {
            if (istart != 0) {
                const double h{sample_ts_[istart] - sample_ts_[istart - 1]};
                const double ratio{(hard_fence.bound_line.minT() - sample_ts_[istart - 1]) / h};
                const std::array<double, 4> proration_coeffs = prorationCoeffs(ratio, h);
                const std::unordered_map<int, double> constrain_vec{
                    {sIndex(istart - 1), proration_coeffs[0]},
                    {sIndex(istart), proration_coeffs[1]},
                    {vIndex(istart - 1), proration_coeffs[2]},
                    {vIndex(istart), proration_coeffs[3]}};
                qp_problem_.addConstrainTerm(
                    constrain_vec, hard_fence.bound_line.ys().front(),
                    std::numeric_limits<double>::max()
                );
            }
            for (std::size_t i{istart}; i < iend; ++i) {
                lower_bound[i] = std::max(hard_fence.bound_line(sample_ts_[i]), lower_bound[i]);
            }
            if (iend != num_samples_) {
                const double h{sample_ts_[iend] - sample_ts_[iend - 1]};
                const double ratio{(hard_fence.bound_line.maxT() - sample_ts_[iend - 1]) / h};
                const std::array<double, 4> proration_coeffs = prorationCoeffs(ratio, h);
                const std::unordered_map<int, double> constrain_vec{
                    {sIndex(iend - 1), proration_coeffs[0]},
                    {sIndex(iend), proration_coeffs[1]},
                    {vIndex(iend - 1), proration_coeffs[2]},
                    {vIndex(iend), proration_coeffs[3]}};
                qp_problem_.addConstrainTerm(
                    constrain_vec, hard_fence.bound_line.ys().back(),
                    std::numeric_limits<double>::max()
                );
            }
        } else {
            if (istart != 0) {
                const double h{sample_ts_[istart] - sample_ts_[istart - 1]};
                const double ratio{(hard_fence.bound_line.minT() - sample_ts_[istart - 1]) / h};
                const std::array<double, 4> proration_coeffs = prorationCoeffs(ratio, h);
                const std::unordered_map<int, double> constrain_vec{
                    {sIndex(istart - 1), proration_coeffs[0]},
                    {sIndex(istart), proration_coeffs[1]},
                    {vIndex(istart - 1), proration_coeffs[2]},
                    {vIndex(istart), proration_coeffs[3]}};
                qp_problem_.addConstrainTerm(
                    constrain_vec, std::numeric_limits<double>::lowest(),
                    hard_fence.bound_line.ys().front()
                );
            }
            for (std::size_t i{istart}; i < iend; ++i) {
                upper_bound[i] = std::min(hard_fence.bound_line(sample_ts_[i]), upper_bound[i]);
            }
            if (iend != num_samples_) {
                const double h{sample_ts_[iend] - sample_ts_[iend - 1]};
                const double ratio{(hard_fence.bound_line.maxT() - sample_ts_[iend - 1]) / h};
                const std::array<double, 4> proration_coeffs = prorationCoeffs(ratio, h);
                const std::unordered_map<int, double> constrain_vec{
                    {sIndex(iend - 1), proration_coeffs[0]},
                    {sIndex(iend), proration_coeffs[1]},
                    {vIndex(iend - 1), proration_coeffs[2]},
                    {vIndex(iend), proration_coeffs[3]}};
                qp_problem_.addConstrainTerm(
                    constrain_vec, std::numeric_limits<double>::lowest(),
                    hard_fence.bound_line.ys().back()
                );
            }
        }
    }

    for (std::size_t i{1}; i < num_samples_ - 1; ++i) {
        qp_problem_.updateConstrainTerm(
            sIndex(i), {{sIndex(i), 1.0}}, lower_bound[i], upper_bound[i]
        );
    }
    return;
}

void RouteLineCubicAccModel::setVelocityRange(double lower_bound, double upper_bound) noexcept {
    for (std::size_t i{1}; i < num_samples_ - 1; ++i) {
        qp_problem_.updateConstrainTerm(vIndex(i), {{vIndex(i), 1.0}}, lower_bound, upper_bound);
    }
    return;
}

void RouteLineCubicAccModel::setAccelRange(double lower_bound, double upper_bound) noexcept {
    qp_problem_.updateConstrainTerm(
        num_samples_ * 2,
        {{sIndex(0), -reciprocal_h2s_[0] * 6.0},
         {sIndex(1), reciprocal_h2s_[0] * 6.0},
         {vIndex(0), -reciprocal_hs_[0] * 4.0},
         {vIndex(1), -reciprocal_hs_[0] * 2.0}},
        std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()
    );
    for (std::size_t i{1}; i < num_samples_ - 1; ++i) {
        qp_problem_.updateConstrainTerm(
            num_samples_ * 2 + i,
            {{sIndex(i), -reciprocal_h2s_[i] * 6.0},
             {sIndex(i + 1), reciprocal_h2s_[i] * 6.0},
             {vIndex(i), -reciprocal_hs_[i] * 4.0},
             {vIndex(i + 1), -reciprocal_hs_[i] * 2.0}},
            lower_bound, upper_bound
        );
    }
    qp_problem_.updateConstrainTerm(
        num_samples_ * 3 - 1,
        {{sIndex(num_samples_ - 2), reciprocal_h2s_[num_samples_ - 2] * 6.0},
         {sIndex(num_samples_ - 1), -reciprocal_h2s_[num_samples_ - 2] * 6.0},
         {vIndex(num_samples_ - 2), reciprocal_hs_[num_samples_ - 2] * 2.0},
         {vIndex(num_samples_ - 1), reciprocal_hs_[num_samples_ - 2] * 4.0}},
        std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()
    );
    return;
}

void RouteLineCubicAccModel::setInitialState(double s0, double v0, double a0) noexcept {
    if (!std::isnan(s0)) {
        qp_problem_.updateConstrainTerm(
            sIndex(0), {{sIndex(0), 1.0}}, s0 - ::boyle::math::kEpsilon,
            s0 + ::boyle::math::kEpsilon
        );
    } else {
        qp_problem_.updateConstrainTerm(
            sIndex(0), {{sIndex(0), 1.0}}, std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::max()
        );
    }
    if (!std::isnan(v0)) {
        qp_problem_.updateConstrainTerm(
            vIndex(0), {{vIndex(0), 1.0}}, v0 - ::boyle::math::kEpsilon,
            v0 + ::boyle::math::kEpsilon
        );
    } else {
        qp_problem_.updateConstrainTerm(
            vIndex(0), {{vIndex(0), 1.0}}, std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::max()
        );
    }
    if (!std::isnan(a0)) {
        qp_problem_.updateConstrainTerm(
            num_samples_ * 3,
            {{sIndex(0), -reciprocal_h2s_[0] * 6.0},
             {sIndex(1), reciprocal_h2s_[0] * 6.0},
             {vIndex(0), -reciprocal_hs_[0] * 4.0},
             {vIndex(1), -reciprocal_hs_[0] * 2.0}},
            a0 - ::boyle::math::kEpsilon, a0 + ::boyle::math::kEpsilon
        );
    } else {
        qp_problem_.updateConstrainTerm(
            num_samples_ * 3,
            {{sIndex(0), -reciprocal_h2s_[0] * 6.0},
             {sIndex(1), reciprocal_h2s_[0] * 6.0},
             {vIndex(0), -reciprocal_hs_[0] * 4.0},
             {vIndex(1), -reciprocal_hs_[0] * 2.0}},
            std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()
        );
    }
    return;
}

void RouteLineCubicAccModel::setFinalState(double sf, double vf, double af) noexcept {
    if (!std::isnan(sf)) {
        qp_problem_.updateConstrainTerm(
            sIndex(num_samples_ - 1), {{sIndex(num_samples_ - 1), 1.0}},
            sf - ::boyle::math::kEpsilon, sf + ::boyle::math::kEpsilon
        );
    } else {
        qp_problem_.updateConstrainTerm(
            sIndex(num_samples_ - 1), {{sIndex(num_samples_ - 1), 1.0}},
            std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()
        );
    }
    if (!std::isnan(vf)) {
        qp_problem_.updateConstrainTerm(
            vIndex(num_samples_ - 1), {{vIndex(num_samples_ - 1), 1.0}},
            vf - ::boyle::math::kEpsilon, vf + ::boyle::math::kEpsilon
        );
    } else {
        qp_problem_.updateConstrainTerm(
            vIndex(num_samples_ - 1), {{vIndex(num_samples_ - 1), 1.0}},
            std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()
        );
    }
    if (!std::isnan(af)) {
        qp_problem_.updateConstrainTerm(
            num_samples_ * 4 - 1,
            {{sIndex(num_samples_ - 2), reciprocal_h2s_[num_samples_ - 2] * 6.0},
             {sIndex(num_samples_ - 1), -reciprocal_h2s_[num_samples_ - 2] * 6.0},
             {vIndex(num_samples_ - 2), reciprocal_hs_[num_samples_ - 2] * 2.0},
             {vIndex(num_samples_ - 1), reciprocal_hs_[num_samples_ - 2] * 4.0}},
            af - ::boyle::math::kEpsilon, af + ::boyle::math::kEpsilon
        );
    } else {
        qp_problem_.updateConstrainTerm(
            num_samples_ * 4 - 1,
            {{sIndex(num_samples_ - 2), reciprocal_h2s_[num_samples_ - 2] * 6.0},
             {sIndex(num_samples_ - 1), -reciprocal_h2s_[num_samples_ - 2] * 6.0},
             {vIndex(num_samples_ - 2), reciprocal_hs_[num_samples_ - 2] * 2.0},
             {vIndex(num_samples_ - 1), reciprocal_hs_[num_samples_ - 2] * 4.0}},
            std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max()
        );
    }
    return;
}

void RouteLineCubicAccModel::setVelocityCost(
    double target_velocity, double velocity_weight
) noexcept {
    for (std::size_t i{0}; i < num_samples_ - 1; ++i) {
        const double factor = velocity_weight * reciprocal_hs_[i] / time_scale_;
        qp_problem_.addQuadCostTerm(sIndex(i), sIndex(i), factor * 6.0 / 5.0);
        qp_problem_.addQuadCostTerm(sIndex(i), sIndex(i + 1), -factor * 12.0 / 5.0);
        qp_problem_.addQuadCostTerm(sIndex(i), vIndex(i), factor * hs_[i] / 5.0);
        qp_problem_.addQuadCostTerm(sIndex(i), vIndex(i + 1), factor * hs_[i] / 5.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), sIndex(i + 1), factor * 6.0 / 5.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), vIndex(i), -factor * hs_[i] / 5.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), vIndex(i + 1), -factor * hs_[i] / 5.0);
        qp_problem_.addQuadCostTerm(vIndex(i), vIndex(i), factor * h2s_[i] * 2.0 / 15.0);
        qp_problem_.addQuadCostTerm(vIndex(i), vIndex(i + 1), -factor * h2s_[i] / 15.0);
        qp_problem_.addQuadCostTerm(vIndex(i + 1), vIndex(i + 1), factor * h2s_[i] * 2.0 / 15.0);
        qp_problem_.addLinCostTerm(sIndex(i), factor * target_velocity * hs_[i] * 2.0);
        qp_problem_.addLinCostTerm(sIndex(i + 1), -factor * target_velocity * hs_[i] * 2.0);
    }
    return;
}

void RouteLineCubicAccModel::setAccelCost(double accel_weight) noexcept {
    for (std::size_t i{0}; i < num_samples_ - 1; ++i) {
        const double factor = accel_weight * reciprocal_h2s_[i] * reciprocal_hs_[i] / time_scale_;
        qp_problem_.addQuadCostTerm(sIndex(i), sIndex(i), factor * 12.0);
        qp_problem_.addQuadCostTerm(sIndex(i), sIndex(i + 1), -factor * 24.0);
        qp_problem_.addQuadCostTerm(sIndex(i), vIndex(i), factor * hs_[i] * 12.0);
        qp_problem_.addQuadCostTerm(sIndex(i), vIndex(i + 1), factor * hs_[i] * 12.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), sIndex(i + 1), factor * 12.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), vIndex(i), -factor * hs_[i] * 12.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), vIndex(i + 1), -factor * hs_[i] * 12.0);
        qp_problem_.addQuadCostTerm(vIndex(i), vIndex(i), factor * h2s_[i] * 4.0);
        qp_problem_.addQuadCostTerm(vIndex(i), vIndex(i + 1), factor * h2s_[i] * 4.0);
        qp_problem_.addQuadCostTerm(vIndex(i + 1), vIndex(i + 1), factor * h2s_[i] * 4.0);
    }
    return;
}

void RouteLineCubicAccModel::setJerkCost(double jerk_weight) noexcept {
    for (std::size_t i{0}; i < num_samples_ - 1; ++i) {
        const double factor =
            jerk_weight * reciprocal_h2s_[i] * reciprocal_h2s_[i] * reciprocal_hs_[i] / time_scale_;
        qp_problem_.addQuadCostTerm(sIndex(i), sIndex(i), factor * 144.0);
        qp_problem_.addQuadCostTerm(sIndex(i), sIndex(i + 1), -factor * 288.0);
        qp_problem_.addQuadCostTerm(sIndex(i), vIndex(i), factor * hs_[i] * 144.0);
        qp_problem_.addQuadCostTerm(sIndex(i), vIndex(i + 1), factor * hs_[i] * 144.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), sIndex(i + 1), factor * 144.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), vIndex(i), -factor * hs_[i] * 144.0);
        qp_problem_.addQuadCostTerm(sIndex(i + 1), vIndex(i + 1), -factor * hs_[i] * 144.0);
        qp_problem_.addQuadCostTerm(vIndex(i), vIndex(i), factor * h2s_[i] * 36.0);
        qp_problem_.addQuadCostTerm(vIndex(i), vIndex(i + 1), factor * h2s_[i] * 72.0);
        qp_problem_.addQuadCostTerm(vIndex(i + 1), vIndex(i + 1), factor * h2s_[i] * 36.0);
    }
    return;
}

RouteLineCubicAccModel::Result RouteLineCubicAccModel::solve() const noexcept {
    ::boyle::math::OsqpSolver osqp_solver{settings_};
    ::boyle::math::OsqpResult osqp_result = osqp_solver.solve(qp_problem_);
    const double a0 =
        -(osqp_result.prim_vars[sIndex(0)] - osqp_result.prim_vars[sIndex(1)]) *
            reciprocal_h2s_[0] * 6.0 -
        (osqp_result.prim_vars[vIndex(0)] * 4.0 + osqp_result.prim_vars[vIndex(1)] * 2.0) *
            reciprocal_hs_[0];
    const double af = -(osqp_result.prim_vars[sIndex(num_samples_ - 1)] -
                        osqp_result.prim_vars[sIndex(num_samples_ - 2)]) *
                          reciprocal_h2s_[num_samples_ - 2] * 6.0 +
                      (osqp_result.prim_vars[vIndex(num_samples_ - 1)] * 4.0 +
                       osqp_result.prim_vars[vIndex(num_samples_ - 2)] * 2.0) *
                          reciprocal_hs_[num_samples_ - 2];
    std::array<Motion1d::BoundaryMode, 2> b0{
        Motion1d::BoundaryMode{2, a0}, Motion1d::BoundaryMode{4, 0.0}};
    std::array<Motion1d::BoundaryMode, 2> bf{
        Motion1d::BoundaryMode{2, af}, Motion1d::BoundaryMode{4, 0.0}};
    return Result{
        .motion =
            Motion1d{
                sample_ts_,
                {osqp_result.prim_vars.cbegin(), osqp_result.prim_vars.cbegin() + num_samples_},
                b0,
                bf},
        .info = osqp_result.info};
}

void RouteLineCubicAccModel::clear() noexcept {
    num_samples_ = 0;
    time_scale_ = 0.0;
    sample_ts_.clear();
    hs_.clear();
    h2s_.clear();
    reciprocal_hs_.clear();
    reciprocal_h2s_.clear();
}

std::size_t RouteLineCubicAccModel::sIndex(std::size_t t_index) const noexcept { return t_index; }

std::size_t RouteLineCubicAccModel::vIndex(std::size_t t_index) const noexcept {
    return num_samples_ + t_index;
}

} // namespace kinetics
} // namespace boyle
