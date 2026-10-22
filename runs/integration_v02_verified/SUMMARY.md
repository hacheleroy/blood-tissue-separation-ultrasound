# Benchmark run

This run is a software check; it does not establish clinical or acoustic performance.

| Method | Status | Velocity NRMSE | Near-wall NRMSE | WSS NRMSE | Time (s) |
|---|---|---:|---:|---:|---:|
| svd_fixed_uncentered | ok | 0.0136 | 0.0247 | 1.11 | 0.0234 |
| svd_energy_uncentered | ok | 0.0136 | 0.0247 | 1.11 | 0.0254 |
| svd_baranger_uncentered | ok | 0.0142 | 0.0286 | 0.479 | 0.0392 |
| svd_long_uncentered | ok | 0.013 | 0.0217 | 1.03 | 0.0487 |
| svd_local_uncentered | ok | 0.0162 | 0.0269 | 0.814 | 0.0451 |
| song_local_reconstruction_uncentered | ok | 0.0436 | 0.0771 | 0.37 | 0.276 |
| svd_temporal_frequency_uncentered | ok | 0.0136 | 0.0247 | 1.11 | 0.0299 |
| evd_fixed_uncentered | ok | 0.0136 | 0.0247 | 1.11 | 0.0254 |
| song_rsvd_uncentered | ok | 0.0136 | 0.0247 | 1.11 | 0.0205 |
| song_rsvd_downsampled_uncentered | ok | 0.0136 | 0.0247 | 0.849 | 0.0183 |
| huang_evd_median_uncentered | ok | 4.3 | 6.79 | 1.76 | 0.0257 |
| huang_evd_ncc_uncentered | ok | 0.0648 | 0.0841 | 1.34 | 0.0292 |
| huang_evd_sigmoid_uncentered | ok | 0.0534 | 0.06 | 1.39 | 0.0288 |
| huang_gaussian_iq_extension_uncentered | ok | 0.0544 | 0.0712 | 1.48 | 0.0741 |
| aer_temporal_alpha0.5_uncentered | ok | 0.0295 | 0.0315 | 0.112 | 0.0531 |
| aer_temporal_alpha0.25_uncentered | ok | 0.365 | 0.284 | 0.53 | 0.0523 |
| aer_spatial_alpha0.5_uncentered | ok | 0.0295 | 0.0315 | 0.112 | 0.0416 |
| aer_spatial_alpha0.25_uncentered | ok | 0.365 | 0.284 | 0.53 | 0.0426 |
| aer_noise_floor_uncentered | ok | 0.287 | 0.335 | 1.55 | 0.0492 |
| axis1_uncentered | ok | 0.016 | 0.0254 | 0.714 | 0.0471 |
| axis1_rigid_uncentered | ok | 0.016 | 0.0254 | 0.714 | 0.114 |
| axis2_uncentered | ok | 0.581 | 0.534 | 0.83 | 0.645 |
| axis2_with_memory_uncentered | ok | 0.622 | 0.556 | 0.789 | 0.397 |
| svd_fixed_centered | ok | 0.0148 | 0.03 | 1.26 | 0.0208 |
| svd_energy_centered | ok | 0.0148 | 0.03 | 1.26 | 0.02 |
| svd_baranger_centered | ok | 0.0157 | 0.0349 | 0.51 | 0.0249 |
| svd_long_centered | ok | 0.0132 | 0.0229 | 1.05 | 0.041 |
| svd_local_centered | ok | 0.0172 | 0.0323 | 1 | 0.0544 |
| song_local_reconstruction_centered | ok | 0.0449 | 0.0827 | 0.384 | 0.253 |
| svd_temporal_frequency_centered | ok | 0.629 | 0.501 | 0.804 | 0.022 |
| evd_fixed_centered | ok | 0.0148 | 0.03 | 1.26 | 0.0213 |
| song_rsvd_centered | ok | 0.0148 | 0.03 | 1.26 | 0.0177 |
| song_rsvd_downsampled_centered | ok | 0.0148 | 0.03 | 0.985 | 0.0204 |
| huang_evd_median_centered | ok | 4.31 | 6.82 | 1.73 | 0.0206 |
| huang_evd_ncc_centered | ok | 0.0477 | 0.0794 | 1.44 | 0.0235 |
| huang_evd_sigmoid_centered | ok | 0.0434 | 0.0595 | 1.49 | 0.0232 |
| huang_gaussian_iq_extension_centered | ok | 0.0571 | 0.0617 | 1.57 | 0.0679 |
| aer_temporal_alpha0.5_centered | ok | 0.0291 | 0.0316 | 0.113 | 0.0533 |
| aer_temporal_alpha0.25_centered | ok | 0.362 | 0.282 | 0.527 | 0.046 |
| aer_spatial_alpha0.5_centered | ok | 0.0291 | 0.0316 | 0.113 | 0.0344 |
| aer_spatial_alpha0.25_centered | ok | 0.362 | 0.282 | 0.527 | 0.0358 |
| aer_noise_floor_centered | ok | 0.293 | 0.336 | 1.57 | 0.0436 |
| axis1_centered | ok | 0.0264 | 0.0597 | 1.5 | 0.0465 |
| axis1_rigid_centered | ok | 0.0264 | 0.0597 | 1.5 | 0.106 |
| axis2_centered | ok | 0.0203 | 0.0393 | 0.393 | 1.21 |
| axis2_with_memory_centered | ok | 0.0238 | 0.0415 | 0.382 | 0.633 |
| ica_circular_complex | ok | 0.695 | 0.663 | 0.845 | 2.2 |
| rpca_pcp | ok | 0.159 | 0.181 | 1.41 | 0.703 |
| butterworth_hp_0.01_prf | ok | 0.714 | 0.562 | 0.811 | 0.0222 |
| butterworth_hp_0.03_prf | ok | 0.0415 | 0.0974 | 0.281 | 0.0206 |
| butterworth_hp_0.08_prf | ok | 0.356 | 1.07 | 2.87 | 0.0182 |
| butterworth_hp_causal_0.03_prf | ok | 0.0771 | 0.0784 | 0.192 | 0.0239 |
| oracle_svd_velocity_centeredFalse | ok | 0.0136 | 0.0247 | 1.11 | 0.117 |
| oracle_svd_velocity_centeredTrue | ok | 0.0148 | 0.03 | 1.26 | 0.118 |
| oracle_aer_velocity | ok | 0.0295 | 0.0315 | 0.112 | 0.181 |
