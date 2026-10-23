# Shared three-axis benchmark v05

40 synthetic cases, 2920 case-setting evaluations, three fits per setting. 24 held-out cases.
Frozen source: `c1e959c03e87c005932c4b35938d3ce3d5a10c7b9b4266f3861bc9a00ad47e93`.

Global operating points pool all training conditions. Means below pool the 24 held-out cases; a mean training target is not a per-case guarantee. Known-motion arms and truth-fitted Axis 2 weights are excluded from selection.

| Family | Motion | History | Training target | Selected setting | Held-out mean IQ error | Mean background | Pass fraction |
|---|---|---:|---:|---|---:|---:|---:|
| axis1 | none | 16 | unconstrained | axis1_none_h16_k1 | 4.228 | 0.02816 | — |
| axis1 | none | 16 | 0.01 | infeasible | — | — | — |
| axis1 | none | 16 | 0.0001 | infeasible | — | — | — |
| axis1 | none | 64 | unconstrained | axis1_none_h64_k4 | 2.931 | 0.01062 | — |
| axis1 | none | 64 | 0.01 | axis1_none_h64_k4 | 2.931 | 0.01062 | 0.625 |
| axis1 | none | 64 | 0.0001 | infeasible | — | — | — |
| axis1 | rigid | 64 | unconstrained | axis1_rigid_h64_k1 | 7.486 | 0.04505 | — |
| axis1 | rigid | 64 | 0.01 | infeasible | — | — | — |
| axis1 | rigid | 64 | 0.0001 | infeasible | — | — | — |
| axis1 | smooth | 64 | unconstrained | axis1_smooth_h64_k4 | 5.758 | 0.02188 | — |
| axis1 | smooth | 64 | 0.01 | infeasible | — | — | — |
| axis1 | smooth | 64 | 0.0001 | infeasible | — | — | — |
| axis2 | none | 16 | unconstrained | axis2_none_h16_r4_p0.1 | 24.26 | 0.3151 | — |
| axis2 | none | 16 | 0.01 | infeasible | — | — | — |
| axis2 | none | 16 | 0.0001 | infeasible | — | — | — |
| axis2 | none | 64 | unconstrained | axis2_none_h64_r4_p0.3 | 22.25 | 0.2559 | — |
| axis2 | none | 64 | 0.01 | infeasible | — | — | — |
| axis2 | none | 64 | 0.0001 | infeasible | — | — | — |
| axis2 | rigid | 64 | unconstrained | axis2_rigid_h64_r4_p0.1 | 13.52 | 0.09964 | — |
| axis2 | rigid | 64 | 0.01 | infeasible | — | — | — |
| axis2 | rigid | 64 | 0.0001 | infeasible | — | — | — |
| axis2 | smooth | 64 | unconstrained | axis2_smooth_h64_r4_p0.3 | 12.04 | 0.07726 | — |
| axis2 | smooth | 64 | 0.01 | infeasible | — | — | — |
| axis2 | smooth | 64 | 0.0001 | infeasible | — | — | — |
| baranger | none | 64 | unconstrained | baranger_none_h64 | 0.6816 | 3.869e-06 | — |
| baranger | none | 64 | 0.01 | baranger_none_h64 | 0.6816 | 3.869e-06 | 1 |
| baranger | none | 64 | 0.0001 | infeasible | — | — | — |
| baranger | smooth | 64 | unconstrained | baranger_smooth_h64 | 1.062 | 0.001829 | — |
| baranger | smooth | 64 | 0.01 | baranger_smooth_h64 | 1.062 | 0.001829 | 0.9167 |
| baranger | smooth | 64 | 0.0001 | baranger_smooth_h64 | 1.062 | 0.001829 | 0.5833 |
| lee_spatial_raw | none | 64 | unconstrained | lee_aer_none_h64_spatial_abs26 | 0.9727 | 4.619e-05 | — |
| lee_spatial_raw | none | 64 | 0.01 | lee_aer_none_h64_spatial_abs26 | 0.9727 | 4.619e-05 | 1 |
| lee_spatial_raw | none | 64 | 0.0001 | lee_aer_none_h64_spatial_abs26 | 0.9727 | 4.619e-05 | 0.7083 |
| lee_spatial_raw | smooth | 64 | unconstrained | lee_aer_smooth_h64_spatial_abs26 | 1.036 | 8.398e-05 | — |
| lee_spatial_raw | smooth | 64 | 0.01 | lee_aer_smooth_h64_spatial_abs26 | 1.036 | 8.398e-05 | 1 |
| lee_spatial_raw | smooth | 64 | 0.0001 | lee_aer_smooth_h64_spatial_abs26 | 1.036 | 8.398e-05 | 0.625 |
| lee_temporal_raw | none | 16 | unconstrained | lee_aer_none_h16_temporal_abs26 | 1.048 | 0.0001778 | — |
| lee_temporal_raw | none | 16 | 0.01 | lee_aer_none_h16_temporal_abs26 | 1.048 | 0.0001778 | 1 |
| lee_temporal_raw | none | 16 | 0.0001 | infeasible | — | — | — |
| lee_temporal_raw | none | 64 | unconstrained | lee_aer_none_h64_temporal_abs1 | 1.041 | 0.0003515 | — |
| lee_temporal_raw | none | 64 | 0.01 | lee_aer_none_h64_temporal_abs1 | 1.041 | 0.0003515 | 1 |
| lee_temporal_raw | none | 64 | 0.0001 | lee_aer_none_h64_temporal_abs1000 | 1.022 | 6.47e-05 | 0.625 |
| lee_temporal_raw | rigid | 64 | unconstrained | lee_aer_rigid_h64_temporal_abs26 | 1.027 | 0.0001069 | — |
| lee_temporal_raw | rigid | 64 | 0.01 | lee_aer_rigid_h64_temporal_abs26 | 1.027 | 0.0001069 | 1 |
| lee_temporal_raw | rigid | 64 | 0.0001 | infeasible | — | — | — |
| lee_temporal_raw | smooth | 64 | unconstrained | lee_aer_smooth_h64_temporal_abs1000 | 1.032 | 5.285e-05 | — |
| lee_temporal_raw | smooth | 64 | 0.01 | lee_aer_smooth_h64_temporal_abs1000 | 1.032 | 5.285e-05 | 1 |
| lee_temporal_raw | smooth | 64 | 0.0001 | lee_aer_smooth_h64_temporal_abs1000 | 1.032 | 5.285e-05 | 0.875 |
| lee_unit_extension | none | 64 | unconstrained | lee_aer_none_h64_temporal_unit | 1.365 | 0.0009462 | — |
| lee_unit_extension | none | 64 | 0.01 | lee_aer_none_h64_temporal_unit | 1.365 | 0.0009462 | 1 |
| lee_unit_extension | none | 64 | 0.0001 | infeasible | — | — | — |
| lee_unit_extension | smooth | 64 | unconstrained | lee_aer_smooth_h64_temporal_unit | 1.615 | 0.001004 | — |
| lee_unit_extension | smooth | 64 | 0.01 | lee_aer_smooth_h64_temporal_unit | 1.615 | 0.001004 | 1 |
| lee_unit_extension | smooth | 64 | 0.0001 | infeasible | — | — | — |
| stap_loaded | none | 16 | unconstrained | stap_none_h16_loaded0.01 | 11.66 | 0.22 | — |
| stap_loaded | none | 16 | 0.01 | infeasible | — | — | — |
| stap_loaded | none | 16 | 0.0001 | infeasible | — | — | — |
| stap_loaded | none | 64 | unconstrained | stap_none_h64_loaded0.1 | 9.221 | 0.1783 | — |
| stap_loaded | none | 64 | 0.01 | infeasible | — | — | — |
| stap_loaded | none | 64 | 0.0001 | infeasible | — | — | — |
| stap_loaded | rigid | 64 | unconstrained | stap_rigid_h64_loaded0.01 | 9.633 | 0.07727 | — |
| stap_loaded | rigid | 64 | 0.01 | infeasible | — | — | — |
| stap_loaded | rigid | 64 | 0.0001 | infeasible | — | — | — |
| stap_loaded | smooth | 64 | unconstrained | stap_smooth_h64_loaded0.1 | 13.48 | 0.1068 | — |
| stap_loaded | smooth | 64 | 0.01 | infeasible | — | — | — |
| stap_loaded | smooth | 64 | 0.0001 | infeasible | — | — | — |
| stap_orthogonal | none | 16 | unconstrained | stap_none_h16_orthogonal | 6.095 | 0.05304 | — |
| stap_orthogonal | none | 16 | 0.01 | infeasible | — | — | — |
| stap_orthogonal | none | 16 | 0.0001 | infeasible | — | — | — |
| stap_orthogonal | none | 64 | unconstrained | stap_none_h64_orthogonal | 6.095 | 0.05304 | — |
| stap_orthogonal | none | 64 | 0.01 | infeasible | — | — | — |
| stap_orthogonal | none | 64 | 0.0001 | infeasible | — | — | — |
| stap_orthogonal | smooth | 64 | unconstrained | stap_smooth_h64_orthogonal | 10.57 | 0.06321 | — |
| stap_orthogonal | smooth | 64 | 0.01 | infeasible | — | — | — |
| stap_orthogonal | smooth | 64 | 0.0001 | infeasible | — | — | — |
| svd | none | 16 | unconstrained | svd_none_h16_kall | 1 | 0 | — |
| svd | none | 16 | 0.01 | svd_none_h16_kall | 1 | 0 | 1 |
| svd | none | 16 | 0.0001 | svd_none_h16_kall | 1 | 0 | 1 |
| svd | none | 64 | unconstrained | svd_none_h64_kall | 1 | 0 | — |
| svd | none | 64 | 0.01 | svd_none_h64_kall | 1 | 0 | 1 |
| svd | none | 64 | 0.0001 | svd_none_h64_kall | 1 | 0 | 1 |
| svd | rigid | 64 | unconstrained | svd_rigid_h64_k1 | 8.353 | 0.06093 | — |
| svd | rigid | 64 | 0.01 | infeasible | — | — | — |
| svd | rigid | 64 | 0.0001 | infeasible | — | — | — |
| svd | smooth | 64 | unconstrained | svd_smooth_h64_kall | 1 | 0 | — |
| svd | smooth | 64 | 0.01 | svd_smooth_h64_kall | 1 | 0 | 1 |
| svd | smooth | 64 | 0.0001 | svd_smooth_h64_kall | 1 | 0 | 1 |
| weighted_svd | none | 16 | unconstrained | weighted_svd_none_h16_k1 | 4.228 | 0.02816 | — |
| weighted_svd | none | 16 | 0.01 | infeasible | — | — | — |
| weighted_svd | none | 16 | 0.0001 | infeasible | — | — | — |
| weighted_svd | none | 64 | unconstrained | weighted_svd_none_h64_k4 | 2.931 | 0.01062 | — |
| weighted_svd | none | 64 | 0.01 | weighted_svd_none_h64_k4 | 2.931 | 0.01062 | 0.625 |
| weighted_svd | none | 64 | 0.0001 | infeasible | — | — | — |
| weighted_svd | smooth | 64 | unconstrained | weighted_svd_smooth_h64_k4 | 5.758 | 0.02188 | — |
| weighted_svd | smooth | 64 | 0.01 | infeasible | — | — | — |
| weighted_svd | smooth | 64 | 0.0001 | infeasible | — | — | — |

Maximum frozen closure error: 9.86e-13.
Maximum Axis 1 / matched weighted-SVD primary-metric discrepancy: 2.84e-11.

The accompanying figures and CSVs are regenerated by `scripts/plot_shared_benchmark.py`. Statistical-IQ software results do not establish experimental or acoustic performance.
