# Benchmark run

This run is a software check; it does not establish clinical or acoustic performance.

| Method | Status | Velocity NRMSE | Near-wall NRMSE | WSS NRMSE | Time (s) |
|---|---|---:|---:|---:|---:|
| svd_fixed_uncentered | ok | 0.0129 | 0.0224 | 1.12 | 0.0192 |
| svd_energy_uncentered | ok | 0.0129 | 0.0224 | 1.12 | 0.0195 |
| svd_ssm_surrogate_uncentered | ok | 0.0129 | 0.0224 | 1.12 | 0.0213 |
| svd_long_uncentered | ok | 0.0126 | 0.0214 | 0.991 | 0.0366 |
| svd_local_uncentered | ok | 0.0154 | 0.0228 | 0.836 | 0.043 |
| aer_temporal_alpha0.5_uncentered | ok | 0.0307 | 0.036 | 0.0942 | 0.0419 |
| aer_temporal_alpha0.25_uncentered | ok | 0.374 | 0.38 | 0.76 | 0.0403 |
| aer_spatial_alpha0.5_uncentered | ok | 0.0307 | 0.036 | 0.0942 | 0.0318 |
| aer_spatial_alpha0.25_uncentered | ok | 0.374 | 0.38 | 0.76 | 0.0378 |
| aer_noise_floor_uncentered | ok | 0.265 | 0.349 | 1.61 | 0.0521 |
| axis1_uncentered | ok | 0.0152 | 0.0215 | 0.763 | 0.0405 |
| axis1_rigid_uncentered | ok | 0.0152 | 0.0215 | 0.763 | 0.0872 |
| axis2_uncentered | ok | 0.563 | 0.683 | 1.16 | 0.55 |
| axis2_with_memory_uncentered | ok | 0.601 | 0.7 | 1.11 | 0.41 |
| svd_fixed_centered | ok | 0.013 | 0.0236 | 1.19 | 0.0224 |
| svd_energy_centered | ok | 0.013 | 0.0236 | 1.19 | 0.021 |
| svd_ssm_surrogate_centered | ok | 0.013 | 0.0236 | 1.19 | 0.0224 |
| svd_long_centered | ok | 0.0126 | 0.0212 | 1.03 | 0.0382 |
| svd_local_centered | ok | 0.0153 | 0.0242 | 0.953 | 0.0413 |
| aer_temporal_alpha0.5_centered | ok | 0.0302 | 0.0351 | 0.095 | 0.0405 |
| aer_temporal_alpha0.25_centered | ok | 0.372 | 0.377 | 0.758 | 0.0379 |
| aer_spatial_alpha0.5_centered | ok | 0.0302 | 0.0351 | 0.095 | 0.0334 |
| aer_spatial_alpha0.25_centered | ok | 0.372 | 0.377 | 0.758 | 0.0466 |
| aer_noise_floor_centered | ok | 0.271 | 0.351 | 1.62 | 0.045 |
| axis1_centered | ok | 0.0305 | 0.0797 | 1.43 | 0.0383 |
| axis1_rigid_centered | ok | 0.0305 | 0.0797 | 1.43 | 0.0948 |
| axis2_centered | ok | 0.015 | 0.0284 | 0.328 | 0.828 |
| axis2_with_memory_centered | ok | 0.0157 | 0.0255 | 0.31 | 0.403 |
| ica_circular_complex | ok | 0.709 | 0.769 | 1.14 | 1.21 |
| rpca_pcp | ok | 0.145 | 0.177 | 1.59 | 0.566 |
| oracle_svd_velocity_centeredFalse | ok | 0.0129 | 0.0224 | 1.12 | 0.0926 |
| oracle_svd_velocity_centeredTrue | ok | 0.013 | 0.0236 | 1.19 | 0.0905 |
| oracle_aer_velocity | ok | 0.0307 | 0.036 | 0.0942 | 0.154 |
