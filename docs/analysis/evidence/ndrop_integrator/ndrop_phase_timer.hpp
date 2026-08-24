#ifndef NDROP_PHASE_TIMER_HPP
#define NDROP_PHASE_TIMER_HPP

#include "ndrop_integrator.hpp"

#include <chrono>

namespace ndrop_test {

class PhaseStopwatch {
public:
  using Member = Real AdvancePhaseTimings::*;

  explicit PhaseStopwatch(AdvancePhaseTimings *timings) : timings_(timings) {
    if (timings_ != nullptr) {
      *timings_ = {};
    }
  }

  PhaseStopwatch(const PhaseStopwatch &) = delete;
  PhaseStopwatch &operator=(const PhaseStopwatch &) = delete;

  ~PhaseStopwatch() { record(); }

  void phase(const Member member) {
    if (timings_ == nullptr) {
      return;
    }
    record();
    member_ = member;
    started_ = Clock::now();
  }

private:
  using Clock = std::chrono::steady_clock;

  void record() {
    if (timings_ == nullptr || member_ == nullptr) {
      return;
    }
    const auto stopped = Clock::now();
    timings_->*member_ +=
        std::chrono::duration<Real>(stopped - started_).count();
    member_ = nullptr;
  }

  AdvancePhaseTimings *timings_ = nullptr;
  Member member_ = nullptr;
  Clock::time_point started_{};
};

} // namespace ndrop_test

#endif
