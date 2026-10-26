# Phase 11: prospective adaptive RF result

26 October 2026. The [adaptive candidate](ESDIVA_PHASE11_DEVELOPMENT.md) and
[ten-scene protocol](ESDIVA_PHASE11_PROSPECTIVE_FREEZE.md) were frozen in local
Git commit `5447696d12f62b7b1022997b8cf2e4b85cddbfe2` **before** RF
generation for newly drawn scatterer clouds 3907 and 3919. Its 87-file
[local freeze manifest](../runs/esdiva_phase11_prospective/local_freeze_manifest.json)
records the source and configuration bytes. An attempt to publish that freeze to
GitHub before acquisition was rejected by automatic approval review. After
explicit authorization, the protocol was [published at `15542f9`](https://github.com/hacheleroy/blood-tissue-separation-ultrasound/commit/15542f953a9d0c5b781fd4abb3de6fe5fc2a85ae)
**after** RF generation; the pre-RF freeze therefore has **local, rather
than remotely timestamped, provenance**. No
prospective score or source data was inspected to revise the candidate. The
[thirty source reports, wide-replay audits and scores](../runs/esdiva_phase11_prospective/fresh/),
[decision](../runs/esdiva_phase11_prospective/decision.json) and
[archive manifest](../runs/esdiva_phase11_prospective/manifest.json) retain
the complete new bank; the large RF and IQ checkpoints are local outside Git.

The **frozen prospective numerical gates pass** in all ten declared scenes.
There are no missing or source-ineligible scenes and no numerical gate
failures. The largest selected-minus-no-motion paired added-blood IQ NRMSE
increase over SVD cut 1 and Axis 1 cut 1 is **+0.01206** (limit +0.020).
The largest selected/no-motion refitted-background tissue-power ratio in a
declared moving scene is **0.89285** (limit 0.900). This tissue margin is
only 0.00715. Both selected and no-motion roundtrip errors are at most
**7.85×10⁻¹⁶** (limit 1×10⁻¹⁰). The unchanged uncorrected slow-flow Axis 2
control reaches at most **0.61167** blood NRMSE and **0.15999** tissue leakage
on the four slow-flow scenes (limits 0.700 and 0.180).

| Cloud | Scene | Largest blood NRMSE increase | Largest tissue ratio | Moving-tissue gate |
|---|---|---:|---:|---|
| 3907 | stationary fast | −0.000003 | 0.99988 | Descriptive |
| 3907 | rigid fast | +0.01206 | 0.79072 | Pass |
| 3907 | rigid slow | −0.04407 | 0.78634 | Pass |
| 3907 | smooth slow | −0.05535 | 0.89285 | Pass |
| 3907 | smooth zero | −0.04940 | 0.89530 | Descriptive |
| 3919 | stationary fast | +0.00028 | 1.17416 | Descriptive |
| 3919 | rigid fast | −0.00313 | 0.71111 | Pass |
| 3919 | rigid slow | −0.04881 | 0.71708 | Pass |
| 3919 | smooth slow | −0.04049 | 0.63389 | Pass |
| 3919 | smooth zero | −0.03025 | 0.63515 | Descriptive |

Each row takes the worse value across the two paired methods. Stationary fast
and smooth zero were declared as controls without a moving-tissue suppression
gate. In particular, one stationary SVD tissue ratio is **1.17416**: the
candidate can amplify tissue in that control, despite passing the declared
moving-scene gate. An increase below +0.020 is a *relative* blood-error
constraint, not an absolute flow-recovery guarantee. The largest selected
Axis 1 blood NRMSE is **0.92655** on a slow-flow scene and **0.96248** on a
zero-flow scene, so preservation of all physical flow is not established.

## Input integrity and branch coverage

The RF source additivity error is at most **1.44×10⁻⁷** and all ten scored
8×16 cores replay exactly from the three source RF views. Every scene has
192 wide and 40 core tissue-support pixels. Maximum evaluator-only
blood/tissue power is **0.00651** in the core and **0.00351** in the wide
support, below 0.100. All ten *declared phase-10 prior-calibration* gates
pass; minimum wide/core validation coherence is **0.92574/0.86522**. The
older acquisition-script angular diagnostic fails in eight scenes, but is
explicitly marked legacy-only in the acquisition report and was not the
predeclared phase-11 eligibility gate. Its result has not been silently
reinterpreted as a pass.

The frozen selector chose **hybrid in 10/10 scenes** and reference in 0/10.
Pre-score magnitude validation gain ranges from approximately −0.001 to
0.573, producing phase fractions 0.5–1.0 and shift scales 0.25–0.75 by the
frozen rule. The adaptive hybrid has independent simulation evidence on two
new particle clouds; the natural reference branch has none from this bank.
The five scenes within a cloud share scatterers, and only 16 frames per scene
are scored. This is acoustic simulation, not experimental phantom data.
Measured three-view scanner acquisitions, structural/flow references, matched
controls and independent motion are still required for physical transfer.

The previous phase-10 stress failure remains a separate historical result.
For audit without the large RF checkpoints, run:

```bash
python scripts/verify_esdiva_phase11_development.py
python scripts/verify_esdiva_phase11_prospective.py
```
