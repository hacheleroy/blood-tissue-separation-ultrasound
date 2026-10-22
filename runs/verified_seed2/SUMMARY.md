# Benchmark run

This run is a software check; it does not establish clinical or acoustic performance.

| Method | Status | Velocity NRMSE | Near-wall NRMSE | WSS NRMSE | Time (s) |
|---|---|---:|---:|---:|---:|
| svd_fixed_uncentered | ok | 0.0152 | 0.0308 | 1.05 | 0.0189 |
| svd_energy_uncentered | ok | 0.0152 | 0.0308 | 1.05 | 0.0179 |
| svd_ssm_surrogate_uncentered | ok | 0.0152 | 0.0308 | 1.05 | 0.0199 |
| svd_long_uncentered | ok | 0.0148 | 0.0285 | 0.976 | 0.0336 |
| svd_local_uncentered | ok | 0.0162 | 0.0326 | 0.727 | 0.0427 |
| aer_temporal_alpha0.5_uncentered | ok | 0.0415 | 0.037 | 0.128 | 0.0398 |
| aer_temporal_alpha0.25_uncentered | ok | 0.376 | 0.329 | 0.423 | 0.0381 |
| aer_spatial_alpha0.5_uncentered | ok | 0.0415 | 0.037 | 0.128 | 0.0284 |
| aer_spatial_alpha0.25_uncentered | ok | 0.376 | 0.329 | 0.423 | 0.028 |
| aer_noise_floor_uncentered | ok | 0.26 | 0.33 | 1.61 | 0.038 |
| axis1_uncentered | ok | 0.0159 | 0.0311 | 0.708 | 0.0317 |
| axis1_rigid_uncentered | ok | 0.0159 | 0.0311 | 0.708 | 0.0683 |
| axis2_uncentered | ok | 0.588 | 0.615 | 0.925 | 0.498 |
| axis2_with_memory_uncentered | ok | 0.619 | 0.638 | 0.991 | 0.316 |
| svd_fixed_centered | ok | 0.0159 | 0.034 | 1.15 | 0.0183 |
| svd_energy_centered | ok | 0.0159 | 0.034 | 1.15 | 0.0179 |
| svd_ssm_surrogate_centered | ok | 0.0159 | 0.034 | 1.15 | 0.0188 |
| svd_long_centered | ok | 0.015 | 0.0292 | 1.01 | 0.0337 |
| svd_local_centered | ok | 0.0168 | 0.0357 | 0.84 | 0.0409 |
| aer_temporal_alpha0.5_centered | ok | 0.041 | 0.0367 | 0.125 | 0.0389 |
| aer_temporal_alpha0.25_centered | ok | 0.374 | 0.327 | 0.422 | 0.038 |
| aer_spatial_alpha0.5_centered | ok | 0.041 | 0.0367 | 0.125 | 0.0294 |
| aer_spatial_alpha0.25_centered | ok | 0.374 | 0.327 | 0.422 | 0.029 |
| aer_noise_floor_centered | ok | 0.265 | 0.336 | 1.63 | 0.0378 |
| axis1_centered | ok | 0.0312 | 0.0797 | 1.19 | 0.0323 |
| axis1_rigid_centered | ok | 0.0312 | 0.0797 | 1.19 | 0.0693 |
| axis2_centered | ok | 0.0189 | 0.0432 | 0.308 | 0.697 |
| axis2_with_memory_centered | ok | 0.0195 | 0.043 | 0.312 | 0.439 |
| ica_circular_complex | ok | 0.71 | 0.718 | 1.1 | 1.36 |
| rpca_pcp | ok | 0.138 | 0.155 | 1.71 | 0.685 |
| oracle_svd_velocity_centeredFalse | ok | 0.0152 | 0.0308 | 1.05 | 0.119 |
| oracle_svd_velocity_centeredTrue | ok | 0.0159 | 0.034 | 1.15 | 0.0936 |
| oracle_aer_velocity | ok | 0.0324 | 0.0503 | 0.203 | 0.166 |
