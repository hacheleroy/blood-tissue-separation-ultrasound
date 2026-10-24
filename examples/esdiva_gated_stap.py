"""Fail closed when an eSDIva adaptive STAP fit lacks support or bounded gain.

The decisions inspect only operator diagnostics fitted from observed IQ and a
spatial tissue training annotation. No isolated source or known motion enters
this gate. Its fallback is the existing orthogonal tissue projector on the
same calibration, views and score frames, and it is explicitly recorded.
"""

from btbench.shared import fit_shared


def fit_gated_loaded_stap(context, spec, min_snapshots_per_dimension=2.,
                          max_operator_gain=1.5):
    if (spec.get("method") != "stap" or not spec.get("adaptive") or
            min_snapshots_per_dimension < 1 or max_operator_gain < 1):
        raise ValueError("Need adaptive STAP and positive support/gain limits")
    candidate = fit_shared(context, spec)
    blocks = candidate.diagnostics["blocks"]
    required = int(min_snapshots_per_dimension*blocks[0]["snapshot_dimension"])
    min_support = min(row["training_snapshots_min"] for row in blocks)
    max_gain = max(row["operator_gain_max"] for row in blocks)
    supported = bool(min_support >= required and max_gain <= max_operator_gain)
    decision = {"selected": "adaptive" if supported else "orthogonal_fallback",
                "required_training_snapshots": required,
                "observed_training_snapshots_min": min_support,
                "max_operator_gain_allowed": max_operator_gain,
                "observed_operator_gain_max": max_gain,
                "support_pass": bool(min_support >= required),
                "gain_pass": bool(max_gain <= max_operator_gain),
                "data_only": True}
    if supported:
        return candidate, decision
    return fit_shared(context, {**spec, "adaptive": False}), decision
