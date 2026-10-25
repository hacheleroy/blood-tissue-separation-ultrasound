# eSDIva phase-5 reserved-bank protocol

25 October 2026. The published
[66-file development freeze](../configs/esdiva_phase5_freeze.json) passed
before seeds **2081** and **2093** were opened. This
[reserved-bank manifest](../configs/esdiva_phase5_heldout_protocol.json)
pins the execution adapter and decision rule before either RF simulation.
The adapter differs from the frozen score runner only in admitting these two
reserved seed IDs. Its decision copy changes only the seed list, the recorded
runner filename, and explanatory outcome wording; all numerical filters,
source gates, thresholds, scene definitions, and the score interval are
identical. The [run guard](../examples/esdiva_phase5_run_heldout.py) verifies
both complete source substitutions and the earlier freeze before RF starts.

Each seed receives four Example 22 conditions on the **same scatterer cloud**:
stationary tissue with fast blood; rigidly moving tissue with fast blood;
smoothly moving tissue with blood at 2% speed; and smoothly moving tissue with
zero blood speed. Pinned eSDIva commit, 10% scatterer fraction, the 0°, −5°,
+5° transmit views, 16 stationary disjoint-time calibration frames, 80
evaluation frames, and the last 16 frames after 64 frames of history remain
fixed. The widened 32-row tissue registration input is reconstructed from
source-resolved RF, while the original eight-row IQ tile remains the score
region. Known displacement and isolated sources are evaluator-only.

The [locked decision](../examples/esdiva_phase5_heldout_decide.py) requires
every scene on **both** reserved seeds to pass original and wider angular
input gates, source identity and contamination, and paired roundtrip NRMSE
≤1×10⁻¹⁰. SVD discard 1 and Axis 1: Motion-compensated virtual long ensemble
must each increase paired added-blood IQ NRMSE by no more than 0.020 relative
to no motion in every scene, and their refitted tissue leakage must be no more
than 0.900 of no motion on rigid fast and smooth slow scenes. The separately
fixed rank-4/cap-0.40 Axis 2: Motion-demodulated multi-lag separation setting
must give slow-flow blood NRMSE ≤0.700 and tissue leakage ≤0.180.
Baranger SSM, Stephen Lee AER, global SVD, and Axis 3: Local ultrasound STAP
remain the same reported controls. No adjustment follows any reserved result.

A pass supports moving to the already specified physical phantom protocol.
A failure remains an explicit negative held-out result and is preserved in
full. Neither outcome by itself demonstrates preservation of all slow blood,
an exact geometric displacement field, or physical phantom validity. Actual
phantom validation requires an acquired phantom RF dataset and independent
motion/flow references.
