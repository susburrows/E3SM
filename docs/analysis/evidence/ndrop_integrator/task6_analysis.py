"""Production-facing synthetic Task 6 analysis entry point."""

import json
from pathlib import Path
from typing import Optional, Tuple

from analyze_screening import (
    GroupKey,
    MatchedErrorStatus,
    compare_matched_error_cost,
)
from ndrop_artifacts import load_validated_analysis_input


def run_task6_analysis(
    identity_manifest_path: Path,
    final_manifest_path: Path,
    output_path: Path,
    *,
    group: GroupKey,
    first_axis: Tuple[str, Optional[int]],
    second_axis: Tuple[str, Optional[int]],
    target_error: float,
    require_retained_contract: bool = True,
) -> Path:
    """Validate paired input, compare primary timing, and write one JSON result."""

    paired = load_validated_analysis_input(
        identity_manifest_path,
        final_manifest_path,
        require_retained_contract=require_retained_contract,
    )
    comparison = compare_matched_error_cost(
        paired.dataset,
        group,
        first_axis,
        second_axis,
        target_error,
    )
    result = {"status": comparison.status.value}
    if comparison.status is MatchedErrorStatus.DETERMINATE:
        if comparison.winner_axis is None:
            raise RuntimeError("determinate comparison is missing its winner")
        method, fe_nsubmix = comparison.winner_axis
        winner = {"method": method}
        if fe_nsubmix is not None:
            winner["fe_nsubmix"] = fe_nsubmix
        result["winner"] = winner

    output = Path(output_path)
    payload = {"comparison": result}
    output.write_bytes(
        (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()
    )
    return output
