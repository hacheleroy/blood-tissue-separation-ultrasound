# Phase-6 development split after the phase-5 reserved failure

25 October 2026. [Phase 5](RESULTS_ESDIVA_PHASE5_HELDOUT.md) failed on the
previously reserved seeds 2081 and 2093. Its numerical source, reports,
thresholds and negative decision remain sealed. Those seeds are now **seen**:
they may explain failure modes but cannot tune or validate a replacement.

The [seed declaration](../configs/esdiva_phase6_split.json) precedes any new
RF simulation: 2309 and 2311 are model-design development clouds, 2333 is a
prospective development check opened only after one candidate is fixed, and
2401 and 2411 form a fresh reserve. Each full cloud has the same four
stationary/rigid/smooth slow/smooth zero scenes, source-resolved RF, disjoint
16-frame tissue calibration, 80 evaluation frames, pinned three views,
32-row registration tile and central eight-row score tile.

## Questions and bounded design search

1. Diagnose core cross-angle coherence, first-half versus second-half
   stability, input noise gain and the 32-row vessel-guard tissue region on
   both new design clouds. Keep the phase-5 core gate as a reported control.
   Compare it with a **wide-supported core diagnostic** that requires the
   original ≥0.700 wide-tissue coherence, ≥128 wide pixels, bounded gain ≤4,
   evaluator-only blood/tissue power ≤0.100, and a first-half-selected core
   ROI of ≥40 pixels with bounded gain ≤4. The core coherence and disjoint-time
   validation are reported, even if that core ROI is below 0.700; no isolated
   source or prescribed motion enters the selector. A wide pass is an input
   support hypothesis, not proof of correct filtering.
2. On the same *new development* clouds, compare the unchanged phase-5 local
   phase field with a bounded grid of phase fractions 0.50, 0.65, 0.80 and
   1.00 and Gaussian widths 2 and 4 pixels. The disjoint spatial validation
   gain remains 0.01. Compare SVD cutoff 1 and **Axis 1: Motion-compensated
   virtual long ensemble** cutoff 1 using paired isolated-blood response and
   refitted tissue leakage. Retain the no-motion, exact roundtrip, Baranger
   SSM, Stephen Lee AER, **Axis 2: Motion-demodulated multi-lag separation**
   and **Axis 3: Local ultrasound STAP** controls. Keep evaluator truth out
   of registration and gate APIs.
3. Do not select on a single attractive scene. A candidate must meet all
   original phase-5 paired bounds on both design clouds: blood NRMSE increase
   ≤0.020 for SVD and Axis 1 in every scene, moving-tissue leakage ratio
   ≤0.900 on rigid fast and smooth slow, paired roundtrip ≤1×10⁻¹⁰, and the
   separately fixed slow Axis 2 blood NRMSE ≤0.700 and leakage ≤0.180.
   Source, original/wide input and blood-contamination checks must be audited
   under the chosen input-support rule. Among complete passes, select the
   largest worst-moving-scene tissue-leakage margin; break ties by smaller
   phase fraction, then larger smoothing width. Report all options and
   failures. If none passes, leave the reserve unopened and document the
   negative development result.

Before seed 2333, publish the selected source, config, score code and decision
rule with hashes. Apply them once to all four prospective scenes. Only a full
prospective pass permits a second complete source/results freeze on `main`
before 2401 or 2411 can be opened. A prospective failure returns to declared
development with a *different*, still-unopened confirmation bank; do not
silently recycle 2333 or the old reserves.

This simulation tier cannot supply independent physical optical motion or
flowmeter truth. A separate phase-5-to-phantom bridge must be preregistered
before any physical acquisition; none is yet available.
