# Benchmark run

This run is a software check; it does not establish clinical or acoustic performance.

| Method | Status | Velocity NRMSE | Near-wall NRMSE | WSS NRMSE | Time (s) |
|---|---|---:|---:|---:|---:|
| svd_fixed_uncentered | ok | 0.0136 | 0.0247 | 1.11 | 0.0206 |
| svd_energy_uncentered | ok | 0.0136 | 0.0247 | 1.11 | 0.0201 |
| svd_ssm_surrogate_uncentered | ok | 0.0136 | 0.0247 | 1.11 | 0.0206 |
| svd_long_uncentered | ok | 0.013 | 0.0217 | 1.03 | 0.0378 |
| svd_local_uncentered | ok | 0.0162 | 0.0269 | 0.814 | 0.0406 |
| aer_temporal_alpha0.5_uncentered | ok | 0.0295 | 0.0315 | 0.112 | 0.042 |
| aer_temporal_alpha0.25_uncentered | ok | 0.365 | 0.284 | 0.53 | 0.0403 |
| aer_spatial_alpha0.5_uncentered | ok | 0.0295 | 0.0315 | 0.112 | 0.0299 |
| aer_spatial_alpha0.25_uncentered | ok | 0.365 | 0.284 | 0.53 | 0.0304 |
| aer_noise_floor_uncentered | ok | 0.287 | 0.335 | 1.55 | 0.0403 |
| axis1_uncentered | ok | 0.016 | 0.0254 | 0.714 | 0.0319 |
| axis1_rigid_uncentered | ok | 0.016 | 0.0254 | 0.714 | 0.0687 |
| axis2_uncentered | ok | 0.581 | 0.534 | 0.83 | 0.499 |
| axis2_with_memory_uncentered | ok | 0.622 | 0.556 | 0.789 | 0.304 |
| svd_fixed_centered | ok | 0.0148 | 0.03 | 1.26 | 0.0196 |
| svd_energy_centered | ok | 0.0148 | 0.03 | 1.26 | 0.0207 |
| svd_ssm_surrogate_centered | ok | 0.0148 | 0.03 | 1.26 | 0.0206 |
| svd_long_centered | ok | 0.0132 | 0.0229 | 1.05 | 0.0355 |
| svd_local_centered | ok | 0.0172 | 0.0323 | 1 | 0.0416 |
| aer_temporal_alpha0.5_centered | ok | 0.0291 | 0.0316 | 0.113 | 0.0405 |
| aer_temporal_alpha0.25_centered | ok | 0.362 | 0.282 | 0.527 | 0.0402 |
| aer_spatial_alpha0.5_centered | ok | 0.0291 | 0.0316 | 0.113 | 0.0322 |
| aer_spatial_alpha0.25_centered | ok | 0.362 | 0.282 | 0.527 | 0.0328 |
| aer_noise_floor_centered | ok | 0.293 | 0.336 | 1.57 | 0.0397 |
| axis1_centered | ok | 0.0264 | 0.0597 | 1.5 | 0.032 |
| axis1_rigid_centered | ok | 0.0264 | 0.0597 | 1.5 | 0.0869 |
| axis2_centered | ok | 0.0203 | 0.0393 | 0.393 | 0.711 |
| axis2_with_memory_centered | ok | 0.0238 | 0.0415 | 0.382 | 0.385 |
| ica_circular_complex | ok | 0.695 | 0.663 | 0.845 | 1.18 |
| rpca_pcp | ok | 0.159 | 0.181 | 1.41 | 0.58 |
| oracle_svd_velocity_centeredFalse | ok | 0.0136 | 0.0247 | 1.11 | 0.0939 |
| oracle_svd_velocity_centeredTrue | ok | 0.0148 | 0.03 | 1.26 | 0.095 |
| oracle_aer_velocity | ok | 0.0295 | 0.0315 | 0.112 | 0.174 |
