# Phase 11: adaptive hybrid development on seen RF

26 October 2026. The [phase-10 independent stress bank](ESDIVA_PHASE10_STRESS_RESULT.md) is now **development data**. Its seeds 3607 and 3619, three moving scenes each, supplied a source-resolved screen of 12 fixed hybrid settings: phase fractions {0.5, 0.6, 0.75, 1.0} crossed with shift scales {0.25, 0.5, 0.75}. The [complete scene-level report](../runs/esdiva_phase11_development/report.json) and [manifest](../runs/esdiva_phase11_development/manifest.json) retain every candidate. No newly declared RF cloud was generated for this screen.

**None** of the twelve fixed settings meets the inherited paired gate on all six scenes. The smallest worst tissue ratio is **0.9181** (limit 0.900), and that setting also exceeds the +0.020 paired blood-error ceiling. The two failed phase-10 scenes favor different strengths. A single globally chosen pair is not supported by this grid.

The [phase-11 candidate](../examples/esdiva_phase11_adaptive.py) instead chooses hybrid strengths once from a data-only diagnostic on the **first 64 mixture frames**, before the scored packet. It measures the median held-pixel validation gain of the existing wide-tile magnitude translation. For gain `g`, `w = clip((g − 0.25)/0.20, 0, 1)`, phase fraction `1 − 0.5w`, and nonperiodic unitary shift scale `0.25 + 0.5w`. A stronger measured translation fit gives more weight to spatial transport; a weaker one gives more to local phase. The same mixture-derived strength is used in the paired tissue-plus-noise evaluator refit. The phase-10 source mask, prior calibration gate, reference-mode selector and transport, 64-frame decision interval, SVD/Axis 1/Axis 2 settings and score frames remain unchanged. Isolated blood/tissue and prescribed simulator motion do not select the strength.

The six seen scenes have pre-score translation validation gains **0.247–0.478**. The adaptive rule produces phase fractions **0.500–1.000** and shift scales **0.250–0.750**. Its largest selected-minus-no-motion paired blood NRMSE increase is **+0.01368** (limit +0.020), largest moving-tissue leakage ratio **0.89420** (limit 0.900), and largest roundtrip error **6.42×10⁻¹⁶** (limit 1×10⁻¹⁰). These are **development** scores on only two correlated particle clouds, with a narrow 0.0058 tissue-ratio margin. They cannot validate generalization. The unchanged uncorrected slow Axis 2 control passed in the prior stress archive.

The new prospective protocol declares seeds **3907 and 3919** and five scenes per seed before their RF generation. If any source, calibration, paired blood, moving-tissue, roundtrip or slow Axis 2 gate fails there, the new candidate fails. Natural reference selection is reported separately; an absence leaves that branch unqualified.

```bash
python scripts/verify_esdiva_phase11_development.py
```

The compact verifier checks the source hashes, complete grid, pre-score rule and arithmetic of the archived development result. Full RF replay requires the local, uncommitted source checkpoints.
