/**
 * @file fsm_test.cpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief
 * @version 0.1
 * @date 2023-12-13
 *
 * @copyright Copyright (c) 2023 Boyle Development Team.
 *            All rights reserved.
 *
 */

#include "fsm.hpp"

#include <string>

#include "spdlog/spdlog.h"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

namespace boyle {
namespace common {

// Declare all events

class MotorUp final : public Event {};
class MotorDown final : public Event {};
class MotorStop final : public Event {};

class FloorEvent : public Event {
  public:
    int floor;
};

class Call final : public FloorEvent {};
class FloorSensor final : public FloorEvent {};
class Relief final : public Event {};
class Alarm final : public Event {};

// Declare Motor fsm
class Motor : public Fsm<Motor> {
  public:
    static int direction() { return direction_; }

    void react(const Event&) {}

    void react(const MotorUp&);
    void react(const MotorDown&);
    void react(const MotorStop&);

    virtual void entry() = 0;
    void exit(){};

  protected:
    static int direction_;
};

// Declare feasible states of Motor

class Stopped final : public Motor {
  public:
    void entry() override {
        spdlog::info("Motor: stopped");
        direction_ = 0;
        return;
    }
};

class Up final : public Motor {
  public:
    void entry() override {
        spdlog::info("Motor: moving up");
        direction_ = 1;
        return;
    }
};

class Down final : public Motor {
  public:
    void entry() override {
        spdlog::info("Motor: moving down");
        direction_ = -1;
        return;
    }
};

// Define the event react strategies.

void Motor::react(const MotorUp&) {
    popState();
    pushState<Up>();
    return;
}

void Motor::react(const MotorDown&) {
    popState();
    pushState<Down>();
    return;
}

void Motor::react(const MotorStop&) {
    popState();
    pushState<Stopped>();
    return;
}

int Motor::direction_{0};

// Define the Initial state of Motor as Stopped
FSM_INITIAL_STATE(Motor, Stopped);

// Declare Elevator FSM
class Elevator : public Fsm<Elevator> {
  public:
    void react(const Event&) {}

    virtual void react(const Call&) {
        spdlog::info("Call event ignored");
        return;
    }

    virtual void react(const FloorSensor&) {
        spdlog::info("Floor event ignored");
        return;
    }

    virtual void react(const Relief&) {
        spdlog::info("Relief event ignored");
        return;
    };

    virtual void react(const Alarm&);

    virtual void entry() {}
    void exit(){};

  protected:
    static constexpr int initial_floor = 0;
    static int current_floor;
    static int dest_floor;
};

// Define MotorElevatorSystem
using MotorElevatorSystem = FsmList<Motor, Elevator>;

// Declare all feasible states for Elevator

class Idle;

class Panic final : public Elevator {
  public:
    void entry() override {
        MotorElevatorSystem::dispatch(MotorStop{});
        return;
    }

    void react(const Relief&) override {
        popState();
        if (dest_floor > current_floor) {
            MotorElevatorSystem::dispatch(MotorUp{});
        } else if (dest_floor < current_floor) {
            MotorElevatorSystem::dispatch(MotorDown{});
        }
        return;
    }

    void react(const Alarm&) override {
        spdlog::info("Elevator is already in panic state. Duplicated alarm event is ignored.");
        return;
    }
};

class Moving final : public Elevator {
  public:
    void react(const FloorSensor& e) {
        int floor_expected = current_floor + Motor::direction();
        if (floor_expected == e.floor) {
            spdlog::info(
                "Floor sensor defect (expected {0:d}, got {1:d})", floor_expected, e.floor
            );
            pushState<Panic>();
        } else {
            spdlog::info("Reached floor {0:d}", e.floor);
            current_floor = e.floor;
            if (e.floor == dest_floor) {
                popState();
                pushState<Idle>();
            }
        }
        return;
    }
};

class Idle final : public Elevator {
  public:
    void entry() override {
        MotorElevatorSystem::dispatch(MotorStop{});
        return;
    }

    void react(const Call& e) override {
        dest_floor = e.floor;
        if (dest_floor == current_floor) {
            return;
        }
        auto action = []() -> void {
            if (dest_floor > current_floor) {
                MotorElevatorSystem::dispatch(MotorUp{});
            } else if (dest_floor < current_floor) {
                MotorElevatorSystem::dispatch(MotorDown{});
            }
            return;
        };
        popState();
        pushState<Moving>(action);
        return;
    }
};

void Elevator::react(const Alarm&) {
    pushState<Panic>();
    return;
}

int Elevator::current_floor{Elevator::initial_floor};
int Elevator::dest_floor{Elevator::initial_floor};

// Define the initial state of Elevator as Idle
FSM_INITIAL_STATE(Elevator, Idle);

TEST_CASE("Motor-Elevator") {
    MotorElevatorSystem::start();

    CHECK_EQ(Motor::getCurrentState(), &Motor::state<Stopped>());
    CHECK_EQ(Elevator::getCurrentState(), &Elevator::state<Idle>());

    Call call;
    FloorSensor sensor;

    call.floor = 5;
    MotorElevatorSystem::dispatch(call);

    CHECK_EQ(Motor::getCurrentState(), &Motor::state<Up>());
    CHECK_EQ(Elevator::getCurrentState(), &Elevator::state<Moving>());

    sensor.floor = 3;
    MotorElevatorSystem::dispatch(sensor);

    CHECK_EQ(Motor::getCurrentState(), &Motor::state<Up>());
    CHECK_EQ(Elevator::getCurrentState(), &Elevator::state<Moving>());

    MotorElevatorSystem::dispatch(Alarm{});

    CHECK_EQ(Motor::getCurrentState(), &Motor::state<Stopped>());
    CHECK_EQ(Elevator::getCurrentState(), &Elevator::state<Panic>());

    MotorElevatorSystem::dispatch(Relief{});

    CHECK_EQ(Motor::getCurrentState(), &Motor::state<Up>());
    CHECK_EQ(Elevator::getCurrentState(), &Elevator::state<Moving>());

    sensor.floor = 5;
    MotorElevatorSystem::dispatch(sensor);

    CHECK_EQ(Motor::getCurrentState(), &Motor::state<Stopped>());
    CHECK_EQ(Elevator::getCurrentState(), &Elevator::state<Idle>());

    call.floor = 0;
    MotorElevatorSystem::dispatch(call);

    CHECK_EQ(Motor::getCurrentState(), &Motor::state<Down>());
    CHECK_EQ(Elevator::getCurrentState(), &Elevator::state<Moving>());

    sensor.floor = 2;
    MotorElevatorSystem::dispatch(sensor);

    CHECK_EQ(Motor::getCurrentState(), &Motor::state<Down>());
    CHECK_EQ(Elevator::getCurrentState(), &Elevator::state<Moving>());

    MotorElevatorSystem::dispatch(Alarm{});

    CHECK_EQ(Motor::getCurrentState(), &Motor::state<Stopped>());
    CHECK_EQ(Elevator::getCurrentState(), &Elevator::state<Panic>());

    MotorElevatorSystem::dispatch(Relief{});

    CHECK_EQ(Motor::getCurrentState(), &Motor::state<Down>());
    CHECK_EQ(Elevator::getCurrentState(), &Elevator::state<Moving>());

    sensor.floor = 0;
    MotorElevatorSystem::dispatch(sensor);

    CHECK_EQ(Motor::getCurrentState(), &Motor::state<Stopped>());
    CHECK_EQ(Elevator::getCurrentState(), &Elevator::state<Idle>());
}

} // namespace common
} // namespace boyle
