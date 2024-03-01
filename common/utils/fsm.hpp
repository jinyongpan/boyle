/**
 * @file fsm.hpp
 * @author Houchen Li (houchen_li@hotmail.com)
 * @brief This FSM is derived from digint's tinyfsm project (https://github.com/digint/tinyfsm).
 * @version 0.1
 * @date 2023-12-12
 *
 * @copyright Copyright (c) 2023 Boyle Development Team.
 *            All rights reserved.
 *
 */

#pragma once

#include <functional>
#include <stack>

#include "common/utils/macros.hpp"

namespace boyle {
namespace common {

namespace detail {

template <typename S>
class StateInstance final {
  public:
    using value_type = S;
    using type = StateInstance<S>;
    static S value;
};

template <typename S>
typename StateInstance<S>::value_type StateInstance<S>::value;

} // namespace detail

class Event {
  public:
    Event() noexcept = default;
    DISABLE_COPY_AND_MOVE(Event);
    virtual ~Event() noexcept = 0;
};

inline Event::~Event() noexcept = default;

template <typename F>
class Fsm {
  public:
    Fsm() noexcept = default;
    DISABLE_COPY_AND_MOVE(Fsm);
    virtual ~Fsm() noexcept = 0;

    template <typename S>
    static constexpr S& state() {
        return detail::StateInstance<S>::value;
    }

    static void initialize();

    static void reset() {}

    static void enter() {
        state_stack_.top()->entry();
        return;
    }

    static void start() {
        initialize();
        enter();
        return;
    }

    template <typename E>
    static void dispatch(const E& event) {
        state_stack_.top()->react(event);
        return;
    }

    static F* getCurrentState() { return state_stack_.top(); }

  protected:
    template <typename S>
    static void pushState() {
        state_stack_.top()->exit();
        state_stack_.push(&detail::StateInstance<S>::value);
        state_stack_.top()->entry();
        return;
    }

    template <typename S>
    static void pushState(std::function<void()> action_function) {
        state_stack_.top()->exit();
        action_function();
        state_stack_.push(&detail::StateInstance<S>::value);
        state_stack_.top()->entry();
        return;
    }

    template <typename S>
    static void pushState(
        std::function<void()> action_function, std::function<bool()> condition_function
    ) {
        if (condition_function()) {
            pushState<S>(action_function);
        }
        return;
    }

    static void popState() {
        if (state_stack_.size() >= 1) {
            state_stack_.top()->exit();
            state_stack_.pop();
            if (!state_stack_.empty()) {
                state_stack_.top()->entry();
            }
        }
        return;
    }

  private:
    static std::stack<F*> state_stack_;
};

template <typename F>
inline Fsm<F>::~Fsm() noexcept = default;

template <typename F>
std::stack<F*> Fsm<F>::state_stack_{};

template <typename... FF>
class FsmList;

template <>
class FsmList<> {
  public:
    FsmList() noexcept = default;
    DISABLE_COPY_AND_MOVE(FsmList);
    virtual ~FsmList() noexcept = default;

    static void initialize() {}

    static void reset() {}

    static void enter() {}

    template <typename E>
    static void dispatch(const E&) {}
};

template <typename F, typename... FF>
class [[nodiscard]] FsmList<F, FF...> {
  public:
    FsmList() noexcept = default;
    DISABLE_COPY_AND_MOVE(FsmList);
    virtual ~FsmList() noexcept = default;

    static void initialize() {
        Fsm<F>::initialize();
        FsmList<FF...>::initialize();
        return;
    }

    static void reset() {}

    static void enter() {
        Fsm<F>::enter();
        FsmList<FF...>::enter();
        return;
    }

    static void start() {
        initialize();
        enter();
        return;
    }

    template <typename E>
    static void dispatch(const E& event) {
        Fsm<F>::template dispatch<E>(event);
        FsmList<FF...>::template dispatch<E>(event);
        return;
    }
};

template <typename... SS>
class StateList;

template <>
class StateList<> {
  public:
    StateList() noexcept = default;
    DISABLE_COPY_AND_MOVE(StateList);
    virtual ~StateList() noexcept = default;

    static void reset() {}
};

template <typename S, typename... SS>
class [[nodiscard]] StateList<S, SS...> {
  public:
    StateList() noexcept = default;
    DISABLE_COPY_AND_MOVE(StateList);
    virtual ~StateList() noexcept = default;

    static void reset() {
        S::reset();
        StateList<SS...>::reset();
        return;
    }
};

} // namespace common
} // namespace boyle

#define FSM_INITIAL_STATE(_FSM, _STATE)                                            \
    template <>                                                                    \
    void Fsm<_FSM>::initialize() {                                                 \
        while (!state_stack_.empty()) {                                            \
            state_stack_.pop();                                                    \
        }                                                                          \
        state_stack_.push(&boyle::common::detail::StateInstance<_STATE>::value); \
        return;                                                                    \
    }
