# Phase 10: independent reference and slow-flow RF stress result

26 October 2026. The [primary six-scene protocol](ESDIVA_PHASE10_STRESS_FREEZE.md) was committed at [`f12dc0d`](https://github.com/hacheleroy/blood-tissue-separation-ultrasound/commit/f12dc0d4195030bb4f14eadcb9588a7a2068190c) before RF generation for previously unseen scatterer clouds **3607** and **3619**. The [source-separated reports, wide replay audits, paired scores and reference counterfactuals](../runs/esdiva_phase10_stress/fresh/), [primary decision](../runs/esdiva_phase10_stress/decision.json), [separate reference branch decision](../runs/esdiva_phase10_stress/ablation_decision.json) and [manifest](../runs/esdiva_phase10_stress/manifest.json) are archived. The large RF checkpoints and IQ arrays remain local, outside Git.

The primary decision is **independent stress input or numerical gate failed**: no scene is missing, all six source and phase-10 calibration gates pass, but four paired method/scene tissue-leakage entries fail the fixed ≤0.900 limit. **Natural reference coverage is 0/6**, so this bank does not independently qualify the reference selector. The original 2401/2411 reserve pass is retained as its own narrower result; it is not overwritten or retrospectively relabeled.

## Protocol integrity and input checks

The RF source, three-view geometry (0°, −5°, +5°), 10% scatterer fraction, prior 16-frame calibration, 80 evaluation frames, 32×16 wide tile, 8×16 scored core and fixed phase-10 numerical candidate are unchanged. Six acquisitions cross each cloud with rigid fast (flow scale 1), rigid slow (0.02), and smooth slow (0.02). The first 64 evaluation frames select the mode; only frames 64–79 score. There is no seed, scene or threshold replacement. RF tissue-plus-blood additivity error is at most **1.44×10⁻⁷**. The scored core replays **exactly** from source RF in every scene. Maximum evaluator-only blood/tissue power is **0.0395** in the core and **0.0104** in the wide tissue support, both below 0.100. Minimum wide/core prior-calibration validation coherence is **0.915/0.800**, above the fixed 0.700/0.500 bounds. Three 3619 acquisitions fail the older, separate acquisition-script angular diagnostic; their **actual phase-10** prior-calibration gates pass. Both selected and no-motion roundtrips stay below **6.24×10⁻¹⁶** (limit 1×10⁻¹⁰).

The `comparison.decision` prose inside the inherited parent settings still describes that historical eight-scene reserve. The stress config's explicit six-scene matrix and paired gate rule, applied by its reviewer, define this separate decision.

The initial copied scorer expected a top-level alias for the **already declared 0.1** source-purity limit and raised a `KeyError` before producing a score. The [schema repair](ESDIVA_PHASE10_STRESS_SCHEMA_REPAIR.md) and its [`c487621` commit](https://github.com/hacheleroy/blood-tissue-separation-ultrasound/commit/c487621b4a6dc3ac37ceafe0eec716d79fc538bb) add only that alias and update the counterfactual's parent hash. It was committed **after RF generation started but before any new scene's mode or filter score was inspected**. Thus the semantic threshold remained predeclared, while the executable config bytes were not a strict pre-RF freeze. That provenance limit applies to this result.

## Paired selected-mode outcomes

Each cell is selected-minus-no-motion paired added-blood IQ NRMSE **Δ / selected-to-no-motion refitted-background tissue leakage ratio**. Every scene naturally selected **hybrid**. The fixed requirements are Δ≤+0.020 and tissue ratio≤0.900 for both SVD cut 1 and Axis 1 cut 1, in every scene.

| Cloud | Scene | SVD Δ / tissue ratio | Axis 1 Δ / tissue ratio |
|---:|---|---:|---:|
| 3607 | Rigid fast | +0.0057 / 0.7217 | +0.0001 / 0.7682 |
| 3607 | Rigid slow | −0.0640 / **0.9357** | −0.0521 / **0.9397** |
| 3607 | Smooth slow | −0.0949 / 0.8680 | −0.0723 / 0.8863 |
| 3619 | Rigid fast | +0.0170 / **0.9064** | +0.0080 / **0.9112** |
| 3619 | Rigid slow | −0.0006 / 0.8427 | −0.0045 / 0.8624 |
| 3619 | Smooth slow | +0.0048 / 0.8342 | −0.0004 / 0.8565 |

The maximum blood-error increase is **+0.01703**, within the +0.020 ceiling. Tissue leakage fails for **3607 rigid slow** and **3619 rigid fast**, for both methods. The additional rigid-slow scene probes a condition absent from the earlier four-scene reserve. On the four slow scenes, the unchanged no-motion Axis 2 rank-4/cap-0.4 control gives blood NRMSE **0.560–0.632** (limit 0.700) and tissue leakage **0.145–0.162** (limit 0.180). Passing this control does not erase the selected motion failures or prove that all flow information survives.

## Reference coverage and separate branch sensitivity

The fixed selector requires pre-score local phase loss >0.28, reference-minus-local loss <0.18, and prior core reference support ≥0.9. It requested and selected **hybrid in all six** new scenes; measured core support ranged **0.725–0.825**. No RF scene independently exercised the reference choice. The [separate counterfactual](ESDIVA_PHASE10_REFERENCE_ABLATION_FREEZE.md) was committed at [`2440fdf`](https://github.com/hacheleroy/blood-tissue-separation-ultrasound/commit/2440fdf369869a8a580de6cc713ff1e648afbf80) after RF generation began but before any new score or mode was inspected. It explicitly **forces** the unchanged reference context on all four rigid scenes. This tests branch execution on the new source-resolved RF without changing the primary mode or decision.

| Cloud | Rigid scene | Prior core reference support | Accepted view-frames / 240 | Forced-reference SVD tissue ratio | Forced-reference Axis 1 tissue ratio |
|---:|---|---:|---:|---:|---:|
| 3607 | Fast | 0.825 | 119 | 0.9246 | 0.9395 |
| 3607 | Slow | 0.775 | 115 | 0.9736 | 0.9755 |
| 3619 | Fast | 0.800 | 118 | 1.0306 | 1.0241 |
| 3619 | Slow | 0.725 | 127 | 1.0510 | 1.0423 |

All four forced-reference pairs preserve the paired blood-error-increase ceiling, but all eight forced tissue ratios exceed 0.900; two scenes are even worse than no-motion (>1). These are deliberately **ineligible mode counterfactuals** and cannot count as a natural reference-mode pass. The ordinary mixture fit and selection see no source truth; the paired evaluator separately refits a tissue-plus-noise background using isolated simulator components to measure leakage. The archived `uses_source_truth_for_fit_or_selection: false` flag in the counterfactual refers to the ordinary mixture/selection path, not that evaluator-only background refit.

## Interpretation

This two-cloud RF stress reveals a narrower operating domain than the passed historical reserve, including a moving-tissue failure with fast flow and one with rigid slow flow. It also leaves natural reference selection unvalidated. There are only 16 scored frames per correlated scene, one 10% scatterer tier, a simplified acoustic phantom and prescribed tissue motion. No independent laboratory blood flow, pump control, scanner timing, or tissue-motion reference was measured here. The prior [physical handoff dry-run](PHANTOM_PHASE10_DRYRUN_RESULT.md) remains an input/packet engineering check. A physical study can still acquire and report independent matched captures, but the current numerical candidate cannot be claimed to pass this broader stress or to preserve **all flow**. Any changed candidate needs development on these now-seen clouds and a newly frozen independent test bank before a renewed generalization claim.

```bash
python scripts/verify_esdiva_phase10_stress.py
```

The verifier checks all six primary and four counterfactual records against their pinned source/config hashes and regenerates both decisions byte for byte, without the large RF files.
