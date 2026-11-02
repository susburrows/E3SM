#ifndef NDROP_REPLAY_SCHEMA_HPP
#define NDROP_REPLAY_SCHEMA_HPP

#include "ndrop_integrator.hpp"

#include <string>

namespace ndrop_test::production_replay_v1 {
struct ReferenceSelfCheck;

// Internal, unadmitted measurement. The Python artifact writer supplies and
// verifies run identity, approval references, hashes, and field ordering.
std::string canonical_real_hex(Real value);
std::string serialize_diagnostic_measurement(const AdvanceResult &result,
                                             const std::string &method);
// Closed result-schema component, not a replay-result admission or gate.
std::string serialize_reference_self_check(const ReferenceSelfCheck &check);

} // namespace ndrop_test::production_replay_v1

#endif
