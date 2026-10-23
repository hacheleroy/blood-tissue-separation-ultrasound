# Frozen follow-up execution

35 cases × 110 settings; 3850 returned records.

Two calibration seeds; three held-out seeds per condition. No evaluation truth selects settings.
Primary comparison uses refitted background response and actual output error. Frozen attribution is secondary.

Maximum relative difference of equivalent mixture outputs: 1.35122e-12.
Maximum normalized difference in frozen blood attribution: 0.996049.

| Condition | Family | Selected setting at calibration target 1e-4 | Held-out median output IQ error | Held-out median background/tissue power | Held-out target passes |
|---|---|---|---:|---:|---:|
| slow_overlap | Short SVD | svd_short_c0_k1 | 0.915253 | 2.28678e-06 | 3/3 |
| slow_overlap | Long SVD | svd_long_c0_k1 | 0.352802 | 2.35096e-06 | 3/3 |
| slow_overlap | Local short SVD | local_svd_short_c0_k1 | 0.915346 | 2.21439e-06 | 3/3 |
| slow_overlap | Axis 1, memory | axis1_memory_c0_k1 | 0.613622 | 2.25574e-06 | 3/3 |
| slow_overlap | Axis 1, rigid + memory | axis1_rigid_memory_c0_k1 | 0.613622 | 2.25574e-06 | 3/3 |
| slow_overlap | Lee temporal ZCA | lee_temporal_w256_relative1e-06 | 0.836965 | 3.08032e-06 | 3/3 |
| slow_overlap | Lee spatial ZCA | lee_spatial_absolute1 | 0.76638 | 4.80467e-06 | 3/3 |
| slow_overlap | Axis 2, original | axis2_original_c1_m1 | 0.680251 | 3.01855e-08 | 3/3 |
| slow_overlap | Axis 2, ridge + cap | axis2_ridge_cap_c1_e0.01 | 0.681656 | 3.01869e-08 | 3/3 |
| accelerating | Short SVD | svd_short_c0_k1 | 0.32636 | 2.28827e-05 | 3/3 |
| accelerating | Long SVD | svd_long_c0_k1 | 0.173673 | 2.35183e-05 | 3/3 |
| accelerating | Local short SVD | local_svd_short_c0_k1 | 0.328516 | 2.21563e-05 | 3/3 |
| accelerating | Axis 1, memory | axis1_memory_c0_k1 | 0.276029 | 2.25695e-05 | 3/3 |
| accelerating | Axis 1, rigid + memory | axis1_rigid_memory_c0_k1 | 0.276029 | 2.25695e-05 | 3/3 |
| accelerating | Lee temporal ZCA | lee_temporal_absolute1 | 0.693751 | 3.20029e-05 | 3/3 |
| accelerating | Lee spatial ZCA | lee_spatial_absolute1 | 0.572384 | 4.80555e-05 | 3/3 |
| accelerating | Axis 2, original | axis2_original_c1_m0 | 0.342137 | 7.21962e-07 | 3/3 |
| accelerating | Axis 2, ridge + cap | axis2_ridge_cap_c1_e0.0001 | 0.415777 | 2.95102e-07 | 3/3 |
| texture_jump | Short SVD | svd_short_c0_k1 | 0.621158 | 2.49929e-05 | 3/3 |
| texture_jump | Long SVD | svd_long_c0_k4 | 0.644834 | 2.40184e-05 | 3/3 |
| texture_jump | Local short SVD | local_svd_short_c0_k1 | 0.623816 | 2.41915e-05 | 3/3 |
| texture_jump | Axis 1, memory | axis1_memory_c1_k1 | 0.6776 | 2.35225e-05 | 3/3 |
| texture_jump | Axis 1, rigid + memory | axis1_rigid_memory_c1_k1 | 0.718307 | 2.20271e-05 | 3/3 |
| texture_jump | Lee temporal ZCA | lee_temporal_absolute1 | 0.765772 | 4.46196e-05 | 3/3 |
| texture_jump | Lee spatial ZCA | lee_spatial_absolute1 | 0.678125 | 7.87437e-05 | 3/3 |
| texture_jump | Axis 2, original | axis2_original_c1_m0 | 0.65135 | 7.43544e-07 | 3/3 |
| texture_jump | Axis 2, ridge + cap | infeasible on calibration | n/a | n/a | 0/0 |
| opposed_flow | Short SVD | svd_short_c0_k1 | 0.37251 | 2.28678e-05 | 3/3 |
| opposed_flow | Long SVD | svd_long_c0_k1 | 0.145375 | 2.35096e-05 | 3/3 |
| opposed_flow | Local short SVD | local_svd_short_c0_k1 | 0.376214 | 2.21439e-05 | 3/3 |
| opposed_flow | Axis 1, memory | axis1_memory_c0_k1 | 0.212896 | 2.25574e-05 | 3/3 |
| opposed_flow | Axis 1, rigid + memory | axis1_rigid_memory_c0_k1 | 0.212896 | 2.25574e-05 | 3/3 |
| opposed_flow | Lee temporal ZCA | lee_temporal_absolute1 | 0.752696 | 3.19942e-05 | 3/3 |
| opposed_flow | Lee spatial ZCA | lee_spatial_absolute1 | 0.644757 | 4.80462e-05 | 3/3 |
| opposed_flow | Axis 2, original | axis2_original_c1_m1 | 0.305433 | 3.01745e-07 | 3/3 |
| opposed_flow | Axis 2, ridge + cap | axis2_ridge_cap_c1_e0.01 | 0.306547 | 3.01755e-07 | 3/3 |
| translated_texture | Short SVD | svd_short_c0_k4 | 0.919824 | 7.92603e-05 | 3/3 |
| translated_texture | Long SVD | svd_long_c0_k96 | 1 | 0 | 3/3 |
| translated_texture | Local short SVD | local_svd_short_c0_k4 | 0.943113 | 6.90565e-05 | 3/3 |
| translated_texture | Axis 1, memory | infeasible on calibration | n/a | n/a | 0/0 |
| translated_texture | Axis 1, rigid + memory | axis1_rigid_memory_c1_k16 | 1.00012 | 4.53565e-06 | 3/3 |
| translated_texture | Lee temporal ZCA | infeasible on calibration | n/a | n/a | 0/0 |
| translated_texture | Lee spatial ZCA | infeasible on calibration | n/a | n/a | 0/0 |
| translated_texture | Axis 2, original | infeasible on calibration | n/a | n/a | 0/0 |
| translated_texture | Axis 2, ridge + cap | infeasible on calibration | n/a | n/a | 0/0 |
| psf_deforming | Short SVD | svd_short_c0_k32 | 1 | 0 | 3/3 |
| psf_deforming | Long SVD | svd_long_c0_k96 | 1 | 0 | 3/3 |
| psf_deforming | Local short SVD | local_svd_short_c1_k16 | 1.00055 | 1.58287e-06 | 3/3 |
| psf_deforming | Axis 1, memory | axis1_memory_c1_k16 | 1.00209 | 5.61478e-06 | 3/3 |
| psf_deforming | Axis 1, rigid + memory | axis1_rigid_memory_c0_k16 | 1.00148 | 4.27223e-06 | 3/3 |
| psf_deforming | Lee temporal ZCA | lee_temporal_w256_relative0.1 | 1.01191 | 6.24566e-05 | 3/3 |
| psf_deforming | Lee spatial ZCA | infeasible on calibration | n/a | n/a | 0/0 |
| psf_deforming | Axis 2, original | infeasible on calibration | n/a | n/a | 0/0 |
| psf_deforming | Axis 2, ridge + cap | infeasible on calibration | n/a | n/a | 0/0 |
| equal_static_failure_control | Short SVD | svd_short_c0_k32 | 1 | 0 | 3/3 |
| equal_static_failure_control | Long SVD | svd_long_c0_k96 | 1 | 0 | 3/3 |
| equal_static_failure_control | Local short SVD | local_svd_short_c1_k16 | 1.00048 | 1.08863e-06 | 3/3 |
| equal_static_failure_control | Axis 1, memory | axis1_memory_c1_k16 | 1.00119 | 2.95927e-06 | 3/3 |
| equal_static_failure_control | Axis 1, rigid + memory | axis1_rigid_memory_c1_k16 | 1.00119 | 2.95927e-06 | 3/3 |
| equal_static_failure_control | Lee temporal ZCA | lee_temporal_absolute1000 | 1.00001 | 1.98356e-08 | 3/3 |
| equal_static_failure_control | Lee spatial ZCA | lee_spatial_absolute1000 | 1.00001 | 1.98361e-08 | 3/3 |
| equal_static_failure_control | Axis 2, original | axis2_original_c1_m1 | 1.00009 | 2.2825e-07 | 3/3 |
| equal_static_failure_control | Axis 2, ridge + cap | axis2_ridge_cap_c1_e0.0001 | 1.00009 | 2.28247e-07 | 3/3 |

Infeasible settings are not interpolated. Passing the background-only target does not prove tissue leakage inside F(X).
Absolute AER units and amplitude changes remain in the scores; no gain is fitted to source truth.
All thresholds, held-out failures, source distortions and conditional operator diagnostics are retained in the archive.
