/**
 * @file route_line_cubic_acc_model_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-18
 *
 * @copyright Copyright (c) 2023 Boyle Development Team.
 *            All rights reserved.
 *
 */

#include "kinetics/models/route_line_cubic_acc_model.h"

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "cxxopts.hpp"
#include "matplot/matplot.h"

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest/doctest.h"

#include "math/utils.hpp"

namespace {

bool plot_graph;

} // namespace

namespace boyle {
namespace kinetics {

TEST_CASE("CostTest") {
    auto s_profile = [](double t) -> double {
        return 14.0 + 5.0 * t + 2.8 * t * t + 0.5 * t * t * t;
    };
    auto v_profile = [](double t) -> double { return 5.0 + 5.6 * t + 1.5 * t * t; };

    constexpr std::size_t num_samples = 21;
    const std::vector<double> sample_ts = ::boyle::math::linspace(0.0, 20.0, num_samples);
    RouteLineCubicAccModel route_line_acc_model{std::move(sample_ts)};

    std::vector<double> state_vec(num_samples * 2);
    state_vec.reserve(num_samples * 3);
    for (std::size_t i{0}; i < num_samples; ++i) {
        state_vec[i] = s_profile(sample_ts[i]);
        state_vec[num_samples + i] = v_profile(sample_ts[i]);
    }

    double cost;

    SUBCASE("VelocityCost") {
        route_line_acc_model.setVelocityCost(5.0, 20.0);
        cost = route_line_acc_model.qp_problem().cost(state_vec.cbegin(), state_vec.cend());
        CHECK_EQ(
            cost,
            doctest::Approx(2195626.66667 * 20.0 / 20.0 - 500.0).epsilon(::boyle::math::kEpsilon)
        );
    }

    SUBCASE("AccelCost") {
        route_line_acc_model.setAccelCost(10.0);
        cost = route_line_acc_model.qp_problem().cost(state_vec.cbegin(), state_vec.cend());
        CHECK_EQ(cost, doctest::Approx(31347.2 * 10.0 / 20.0).epsilon(math::kEpsilon));
    }

    SUBCASE("JerkCost") {
        route_line_acc_model.setJerkCost(100.0);
        cost = route_line_acc_model.qp_problem().cost(state_vec.cbegin(), state_vec.cend());
        CHECK_EQ(cost, doctest::Approx(180.0 * 100.0 / 20.0).epsilon(1E-6));
    }

    SUBCASE("SoftFence") {
        ::boyle::math::PiecewiseLinearFunction1d bound_line =
            ::boyle::math::PiecewiseLinearFunction1d{{0.0, 20.0}, {-946.0, 3524.0}};
        std::vector<SoftFence1d> soft_fences{
            SoftFence1d(0, ::boyle::common::Actio::PUSHING, bound_line, 100.0, 10.0)};
        for (std::size_t i{0}; i < num_samples; ++i) {
            state_vec.push_back(std::max(223.5 * sample_ts[i] - 946.0 - state_vec[i], 0.0));
        }
        route_line_acc_model.setSoftFences(soft_fences);
        cost = route_line_acc_model.qp_problem().cost(state_vec.cbegin(), state_vec.cend());
        CHECK_EQ(
            cost,
            doctest::Approx(2966.67 * 100.0 / 20.0 + 1059109.523809519 * 10.0 / 20.0).epsilon(1E-3)
        );
    }
}

TEST_CASE("TrivialScene") {
    constexpr std::size_t num_samples = 21;
    const std::vector<double> sample_ts = ::boyle::math::linspace(0.0, 20.0, num_samples);
    RouteLineCubicAccModel route_line_acc_model{std::move(sample_ts)};
    route_line_acc_model.setInitialState(0.0, 0.0);
    route_line_acc_model.setFinalState(100.0, 0.0);
    route_line_acc_model.setVelocityCost(5.0, 20.0);
    route_line_acc_model.setAccelCost(10.0);
    route_line_acc_model.setJerkCost(100.0);
    route_line_acc_model.setVelocityRange(0.0, 20.0);
    route_line_acc_model.setAccelRange(-5.0, 2.0);
    RouteLineCubicAccModel::Result result = route_line_acc_model.solve();
    const Motion1d& motion = result.motion;

    CHECK_EQ(motion.s(0.0), doctest::Approx(0.0).epsilon(1E-2));
    CHECK_EQ(motion.s(10.0), doctest::Approx(50.0).epsilon(1E1));
    CHECK_EQ(motion.s(20.0), doctest::Approx(100.0).epsilon(1E-2));

    if (plot_graph) {
        using namespace matplot;
        std::vector<double> plot_ts = ::boyle::math::linspace(0.0, 20.0, num_samples);
        std::vector<double> plot_ss;
        plot_ss.reserve(num_samples);
        std::vector<double> plot_vs;
        plot_vs.reserve(num_samples);
        std::vector<double> plot_as;
        plot_as.reserve(num_samples);
        std::vector<double> plot_js;
        plot_js.reserve(num_samples);
        std::vector<double> plot_snaps;
        plot_snaps.reserve(num_samples);
        for (double t : plot_ts) {
            plot_ss.push_back(motion.s(t));
            plot_vs.push_back(motion.velocity(t));
            plot_as.push_back(motion.accel(t));
            plot_js.push_back(motion.jerk(t));
            plot_snaps.push_back(motion.snap(t));
        }

        figure_handle fig = figure();
        fig->size(1600, 1000);
        std::vector<axes_handle> ax;
        for (std::size_t i = 0; i < 5; ++i) {
            ax.emplace_back(subplot(2, 3, i));
            ax[i]->grid(on);
        }

        ax[0]->plot(plot_ts, plot_ss, ".");
        ax[1]->plot(plot_ts, plot_vs, ".");
        ax[2]->plot(plot_ts, plot_as, ".");
        ax[3]->plot(plot_ts, plot_js, ".");
        ax[4]->plot(plot_ts, plot_snaps, ".");
        fig->show();
    }
}

TEST_CASE("FenceScene") {
    constexpr std::size_t num_samples = 21;
    const std::vector<double> sample_ts = ::boyle::math::linspace(0.0, 20.0, num_samples);
    RouteLineCubicAccModel route_line_acc_model{std::move(sample_ts)};
    route_line_acc_model.setInitialState(0.0, 0.0);
    route_line_acc_model.setFinalState(100.0, 0.0);
    route_line_acc_model.setVelocityCost(5.0, 20.0);
    route_line_acc_model.setAccelCost(10.0);
    route_line_acc_model.setJerkCost(100.0);
    route_line_acc_model.setVelocityRange(0.0, 20.0);
    route_line_acc_model.setAccelRange(-5.0, 2.0);

    ::boyle::math::PiecewiseLinearFunction1d bound_line_1{{8.1, 8.9}, {60.0, 60.0}};
    ::boyle::math::PiecewiseLinearFunction1d bound_line_2{{3.2, 3.8}, {20.0, 20.0}};

    Motion1d motion;

    SUBCASE("SoftFences") {
        std::vector<SoftFence1d> soft_fences{
            SoftFence1d(0, ::boyle::common::Actio::BLOCKING, bound_line_1, 1000.0, 10.0),
            SoftFence1d(1, ::boyle::common::Actio::PUSHING, bound_line_2, 1000.0, 10.0)};

        route_line_acc_model.setSoftFences(soft_fences);
        motion = route_line_acc_model.solve().motion;
    }

    SUBCASE("HardFences") {
        std::vector<HardFence1d> hard_fences{
            HardFence1d(0, ::boyle::common::Actio::BLOCKING, bound_line_1),
            HardFence1d(1, ::boyle::common::Actio::PUSHING, bound_line_2)};

        route_line_acc_model.setHardFences(hard_fences);
        motion = route_line_acc_model.solve().motion;

        CHECK_EQ(motion.s(bound_line_1.maxT()), doctest::Approx(60.0).epsilon(1E-1));
        CHECK_EQ(motion.s(bound_line_2.minT()), doctest::Approx(20.0).epsilon(1E-2));
    }

    if (plot_graph) {
        using namespace matplot;
        std::vector<double> plot_ts = ::boyle::math::linspace(0.0, 20.0, num_samples);
        std::vector<double> plot_ss;
        plot_ss.reserve(num_samples);
        std::vector<double> plot_vs;
        plot_vs.reserve(num_samples);
        std::vector<double> plot_as;
        plot_as.reserve(num_samples);
        std::vector<double> plot_js;
        plot_js.reserve(num_samples);
        std::vector<double> plot_snaps;
        plot_snaps.reserve(num_samples);
        for (double t : plot_ts) {
            plot_ss.push_back(motion.s(t));
            plot_vs.push_back(motion.velocity(t));
            plot_as.push_back(motion.accel(t));
            plot_js.push_back(motion.jerk(t));
            plot_snaps.push_back(motion.snap(t));
        }

        figure_handle fig = figure();
        fig->size(1600, 1000);
        std::vector<axes_handle> ax;
        for (size_t i = 0; i < 5; ++i) {
            ax.emplace_back(subplot(2, 3, i));
            ax[i]->grid(on);
            ax[i]->hold(on);
        }

        ax[0]->plot(plot_ts, plot_ss, ".");
        ax[0]->plot(
            {bound_line_1.minT(), bound_line_1.maxT()}, {bound_line_1.minY(), bound_line_1.maxY()},
            "r-"
        );
        ax[0]->plot({bound_line_1.minT(), bound_line_1.minT()}, {bound_line_1.minY(), 100.0}, "r-");
        ax[0]->plot({bound_line_1.maxT(), bound_line_1.maxT()}, {bound_line_1.minY(), 100.0}, "r-");
        ax[0]->plot(
            {bound_line_2.minT(), bound_line_2.maxT()}, {bound_line_2.minY(), bound_line_2.maxY()},
            "r-"
        );
        ax[0]->plot({bound_line_2.minT(), bound_line_2.minT()}, {0.0, bound_line_2.minY()}, "r-");
        ax[0]->plot({bound_line_2.maxT(), bound_line_2.maxT()}, {0.0, bound_line_2.maxY()}, "r-");
        ax[1]->plot(plot_ts, plot_vs, ".");
        ax[2]->plot(plot_ts, plot_as, ".");
        ax[3]->plot(plot_ts, plot_js, ".");
        ax[4]->plot(plot_ts, plot_snaps, ".");
        show();
    }
}

TEST_CASE("HeadwayScene") {
    constexpr std::size_t num_samples = 51;
    const std::vector<double> sample_ts = ::boyle::math::linspace(0.0, 10.0, num_samples);
    RouteLineCubicAccModel route_line_acc_model{std::move(sample_ts)};
    route_line_acc_model.setInitialState(0.0, 20.0, 0.0);
    route_line_acc_model.setVelocityCost(20.0, 20.0);
    route_line_acc_model.setAccelCost(10.0);
    route_line_acc_model.setJerkCost(100.0);
    route_line_acc_model.setVelocityRange(0.0, 20.0);
    route_line_acc_model.setAccelRange(-5.0, 2.0);

    Motion1d motion;
    ::boyle::math::PiecewiseLinearFunction1d bound_line;

    SUBCASE("Yield") {
        bound_line = ::boyle::math::PiecewiseLinearFunction1d{{0.0, 10.0}, {20.0, 120.0}};
        std::vector<SoftFence1d> soft_fences{
            SoftFence1d(0, ::boyle::common::Actio::BLOCKING, bound_line, 1000.0, 10.0)};
        route_line_acc_model.setSoftFences(soft_fences);
        motion = route_line_acc_model.solve().motion;
    }

    SUBCASE("Overtake") {
        bound_line = ::boyle::math::PiecewiseLinearFunction1d{{0.0, 20.0}, {10.0, 210.0}};
        std::vector<SoftFence1d> soft_fences{
            SoftFence1d(0, ::boyle::common::Actio::PUSHING, bound_line, 1000.0, 10.0)};
        route_line_acc_model.setSoftFences(soft_fences);
        motion = route_line_acc_model.solve().motion;
    }

    if (plot_graph) {
        using namespace matplot;
        std::vector<double> plot_ts = ::boyle::math::linspace(0.0, 10.0, num_samples);
        std::vector<double> plot_ss;
        plot_ss.reserve(num_samples);
        std::vector<double> plot_vs;
        plot_vs.reserve(num_samples);
        std::vector<double> plot_as;
        plot_as.reserve(num_samples);
        std::vector<double> plot_js;
        plot_js.reserve(num_samples);
        std::vector<double> plot_snaps;
        plot_snaps.reserve(num_samples);
        for (double t : plot_ts) {
            plot_ss.push_back(motion.s(t));
            plot_vs.push_back(motion.velocity(t));
            plot_as.push_back(motion.accel(t));
            plot_js.push_back(motion.jerk(t));
            plot_snaps.push_back(motion.snap(t));
        }

        figure_handle fig = figure();
        fig->size(1600, 1000);
        std::vector<axes_handle> ax;
        for (std::size_t i = 0; i < 5; ++i) {
            ax.emplace_back(subplot(2, 3, i));
            ax[i]->grid(on);
        }
        ax[0]->plot({bound_line.ts()}, {bound_line.ys()}, "r-");
        ax[0]->plot(plot_ts, plot_ss, "-")->line_width(2.0);
        ax[1]->plot(plot_ts, plot_vs, "-")->line_width(2.0);
        ax[2]->plot(plot_ts, plot_as, "-")->line_width(2.0);
        ax[3]->plot(plot_ts, plot_js, "-")->line_width(2.0);
        ax[4]->plot(plot_ts, plot_snaps, "-")->line_width(2.0);
        show();
    }
}

} // namespace kinetics
} // namespace boyle

int main(int argc, const char* argv[]) {
    cxxopts::Options options(
        "piecewise_linear_function_test", "unit test of PiecewiseLinearFunction class"
    );
    options.add_options()(
        "plot-graph", "plot test graph", cxxopts::value<bool>()->default_value("false")
    );
    cxxopts::ParseResult result = options.parse(argc, argv);
    plot_graph = result["plot-graph"].as<bool>();
    doctest::Context context(argc, argv);
    return context.run();
}
