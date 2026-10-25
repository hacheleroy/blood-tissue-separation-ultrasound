# Current research status — 25 October 2026

The [phase-5 comparison on reserved eSDIva seeds 2081 and 2093](docs/RESULTS_ESDIVA_PHASE5_HELDOUT.md)
has completed and **failed its frozen all-scene gate**. The eight compact
source audits and scores, [locked decision](runs/esdiva_phase5_heldout/decision.json),
and [diagnostic figure](docs/figures/esdiva_phase5_heldout_gates.png) are published.
The reserved seeds are now opened and must not be used for retuning.

The subsequent [phase-6 design study](docs/ESDIVA_PHASE6_DESIGN_RESULT.md)
ran twelve declared calibration and motion trials on fresh development seeds
2309 and 2311, with eight source-resolved scenes and an auditable
[decision](runs/esdiva_phase6_design/decision.json). **No candidate passed**
the unchanged paired blood and moving-tissue gate. The small core calibration
failed on seed 2309 but passed on seed 2311; wider tissue support allowed
scoring on both, yet rigid-motion benefit was not robust across seeds.
Prospective development seed 2333 and fresh reserve 2401/2411 remain unopened.
The [phase-5-to-phantom geometry bridge](docs/PHANTOM_PHASE5_BRIDGE.md) is
preregistered conditionally; its pre-capture card requires a future successful
simulation reserve decision and actual measured hardware geometry.

The subsequent [phase-7 boundary-aware complex reference transport
study](docs/ESDIVA_PHASE7_DESIGN_RESULT.md) froze its source and decision
before generating new RF on development seeds 2503 and 2521. Its
[complete decision](runs/esdiva_phase7_design/decision.json) is **no qualified
candidate**: all four 2503 scenes miss the frozen wide-coherence input gate;
2521 is eligible but rigid motion raises paired blood error and more than
doubles tissue leakage, while slow motion more than doubles leakage. The
source-resolved audits and all sixteen old/new score records are retained.
No phase-7 candidate is released to prospective 2333 or the 2401/2411 reserve.

The root README is included in the earlier 66-file development freeze and
retains its pre-held-out wording as a historical record. This file records the
later outcome without changing that freeze. Physical phantom validation remains
unperformed; the [frozen physical protocol](docs/PHANTOM_PROTOCOL_V07.md)
requires measured acquisitions and independent motion and flow references.
