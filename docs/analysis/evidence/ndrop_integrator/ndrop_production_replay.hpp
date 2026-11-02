#ifndef NDROP_PRODUCTION_REPLAY_HPP
#define NDROP_PRODUCTION_REPLAY_HPP

#include "ndrop_integrator.hpp"

#include <optional>
#include <string>

namespace ndrop_test::production_replay_v1 {

// Offline reference classification only. This interface does not admit a
// production capture or certify a replay case.
struct ReferenceSelfCheck {
  std::string applicability = "required";
  std::string attempt_status = "not_attempted";
  std::string outcome = "not_attempted";
  std::optional<std::string> failure_stage;
  std::optional<std::string> failure_reason;
  Real primary_rtol = 1e-12;
  Real tight_rtol = 1e-13;
  bool primary_completed = false;
  bool tight_completed = false;
  std::optional<Real> maximum_normalized_disagreement;
  Real threshold = 2e-11;
  std::optional<bool> passed;
  std::string message;
};

ReferenceSelfCheck classify_reference_self_check(
    bool candidate_ok, bool replayable, bool stress_sentinel,
    const State &initial, const OperatorData &data, Real dt,
    const std::optional<ReferenceResult> &primary,
    const std::optional<ReferenceResult> &tight);

ReferenceSelfCheck run_reference_self_check(bool candidate_ok,
                                            bool replayable,
                                            bool stress_sentinel,
                                            const State &initial,
                                            const OperatorData &data, Real dt);

} // namespace ndrop_test::production_replay_v1

#endif
