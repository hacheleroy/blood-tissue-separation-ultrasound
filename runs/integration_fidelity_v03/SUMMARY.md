# Benchmark run

Component metrics use the mixture-fitted operator, with no fitted gain or phase correction.
Experimental anatomy and processed velocity maps do not supply independent component truth.

| Method | Status | Blood IQ distortion | Tissue power retention | Blood power retention | Spectral TV | Noise gain | Time (s) |
|---|---|---:|---:|---:|---:|---:|---:|
| svd_fixed_uncentered | ok | 0.124 | 6.29e-07 | 0.985 | 0.014 | 0.985 | 0.04 |
| svd_energy_uncentered | ok | 0.124 | 6.29e-07 | 0.985 | 0.014 | 0.985 | 0.0387 |
| svd_baranger_uncentered | ok | 0.141 | 5.69e-07 | 0.98 | 0.016 | 0.565 | 0.0485 |
| svd_long_uncentered | ok | 0.0678 | 6.09e-07 | 0.995 | 0.00377 | 0.888 | 0.0596 |
| svd_local_uncentered | ok | 0.137 | 2.57e-05 | 0.98 | 0.0146 | 0.754 | 0.0775 |
| song_local_reconstruction_uncentered | ok | 0.683 | 9.36e-06 | 0.366 | 0.278 | 0.126 | 0.397 |
| svd_temporal_frequency_uncentered | ok | 0.124 | 6.29e-07 | 0.985 | 0.014 | 0.985 | 0.0359 |
| evd_fixed_uncentered | ok | 0.124 | 6.29e-07 | 0.985 | 0.014 | 0.985 | 0.0316 |
| song_rsvd_uncentered | ok | 0.124 | 6.29e-07 | 0.985 | 0.014 | 0.985 | 0.0253 |
| song_rsvd_downsampled_uncentered | ok | 0.124 | 1.38e-06 | 0.985 | 0.014 | 0.985 | 0.0242 |
| huang_evd_median_uncentered | ok | 0.996 | 7.74e-08 | 0.00727 | 0.96 | 0.577 | 0.032 |
| huang_evd_ncc_uncentered | ok | 0.842 | 2.64e-07 | 0.291 | 0.459 | 0.922 | 0.0343 |
| huang_evd_sigmoid_uncentered | ok | 0.768 | 4.02e-07 | 0.299 | 0.415 | 0.922 | 0.0348 |
| huang_gaussian_iq_extension_uncentered | ok | 0.835 | 4.37e-07 | 0.275 | 0.492 | 0.92 | 0.0676 |
| aer_temporal_alpha0.5_uncentered | ok | 0.909 | 8.84e-06 | 0.00823 | 0.0158 | 0.00877 | 0.0492 |
| aer_temporal_alpha0.25_uncentered | ok | 0.699 | 0.00297 | 0.0906 | 0.00856 | 0.0935 | 0.0495 |
| aer_spatial_alpha0.5_uncentered | ok | 0.895 | 1.18e-05 | 0.011 | 0.0149 | 0.0117 | 0.0356 |
| aer_spatial_alpha0.25_uncentered | ok | 0.676 | 0.00343 | 0.105 | 0.00747 | 0.108 | 0.0359 |
| aer_noise_floor_uncentered | ok | 0.974 | 1.91e-05 | 1.05 | 0.668 | 53.9 | 0.048 |
| axis1_uncentered | ok | 0.145 | 7.52e-06 | 0.96 | 0.0106 | 0.96 | 0.0461 |
| axis1_rigid_uncentered | ok | 0.145 | 7.52e-06 | 0.96 | 0.0106 | 0.96 | 0.26 |
| axis2_uncentered | ok | 1.36 | 0.246 | 2.59 | 0.0671 | 0.631 | 0.581 |
| axis2_with_memory_uncentered | ok | 0.986 | 0.247 | 1.48 | 0.108 | 0.347 | 0.65 |
| svd_fixed_centered | ok | 0.156 | 6.22e-07 | 0.976 | 0.0206 | 0.969 | 0.0386 |
| svd_energy_centered | ok | 0.156 | 6.22e-07 | 0.976 | 0.0206 | 0.969 | 0.0241 |
| svd_baranger_centered | ok | 0.172 | 5.59e-07 | 0.971 | 0.0227 | 0.52 | 0.0402 |
| svd_long_centered | ok | 0.0741 | 6.08e-07 | 0.995 | 0.00458 | 0.887 | 0.0493 |
| svd_local_centered | ok | 0.166 | 2.56e-05 | 0.971 | 0.021 | 0.748 | 0.0815 |
| song_local_reconstruction_centered | ok | 0.689 | 9.39e-06 | 0.361 | 0.282 | 0.119 | 0.424 |
| svd_temporal_frequency_centered | ok | 0.143 | 0.141 | 0.98 | 0.017 | 0.973 | 0.0526 |
| evd_fixed_centered | ok | 0.156 | 6.22e-07 | 0.976 | 0.0206 | 0.969 | 0.0575 |
| song_rsvd_centered | ok | 0.156 | 6.22e-07 | 0.976 | 0.0206 | 0.969 | 0.0192 |
| song_rsvd_downsampled_centered | ok | 0.156 | 1.35e-06 | 0.976 | 0.0206 | 0.969 | 0.0245 |
| huang_evd_median_centered | ok | 0.997 | 7.29e-08 | 0.00681 | 0.963 | 0.554 | 0.0328 |
| huang_evd_ncc_centered | ok | 0.846 | 2.55e-07 | 0.283 | 0.464 | 0.906 | 0.035 |
| huang_evd_sigmoid_centered | ok | 0.773 | 3.72e-07 | 0.291 | 0.417 | 0.906 | 0.0343 |
| huang_gaussian_iq_extension_centered | ok | 0.84 | 4.27e-07 | 0.268 | 0.492 | 0.904 | 0.0856 |
| aer_temporal_alpha0.5_centered | ok | 0.909 | 8.84e-06 | 0.00831 | 0.0164 | 0.00884 | 0.0547 |
| aer_temporal_alpha0.25_centered | ok | 0.699 | 0.00296 | 0.091 | 0.0092 | 0.0937 | 0.0621 |
| aer_spatial_alpha0.5_centered | ok | 0.895 | 1.18e-05 | 0.0111 | 0.0156 | 0.0117 | 0.0471 |
| aer_spatial_alpha0.25_centered | ok | 0.676 | 0.00341 | 0.105 | 0.00803 | 0.108 | 0.0366 |
| aer_noise_floor_centered | ok | 0.978 | 1.92e-05 | 1.06 | 0.67 | 54.6 | 0.0626 |
| axis1_centered | ok | 0.508 | 6.74e-08 | 0.649 | 0.138 | 0.961 | 0.0789 |
| axis1_rigid_centered | ok | 0.508 | 6.74e-08 | 0.649 | 0.138 | 0.961 | 0.308 |
| axis2_centered | ok | 0.268 | 6.1e-10 | 0.859 | 0.0655 | 0.109 | 1.08 |
| axis2_with_memory_centered | ok | 0.376 | 6.29e-10 | 0.766 | 0.091 | 0.111 | 0.978 |
| ica_circular_complex | ok | 0.574 | 6.65e-08 | 0.56 | 0.209 | 0.922 | 2.1 |
| rpca_pcp | ok | — | — | — | — | — | 0.9 |
| butterworth_hp_0.01_prf | ok | 0.278 | 0.284 | 1.05 | 0.0444 | 1.03 | 0.0397 |
| butterworth_hp_0.03_prf | ok | 0.279 | 0.000168 | 0.902 | 0.0893 | 0.942 | 0.0459 |
| butterworth_hp_0.08_prf | ok | 0.459 | 1.12e-06 | 0.668 | 0.253 | 0.808 | 0.0353 |
| butterworth_hp_causal_0.03_prf | ok | 0.741 | 0.00319 | 0.917 | 0.0715 | 0.937 | 0.044 |
| oracle_svd_preservation_centeredFalse | ok | 0.124 | 6.29e-07 | 0.985 | 0.014 | 0.985 | 0.387 |
| oracle_svd_preservation_centeredTrue | ok | 0.156 | 6.22e-07 | 0.976 | 0.0206 | 0.969 | 0.439 |
| oracle_aer_iq | ok | 0.629 | 9.03e-06 | 0.181 | 0.255 | 0.781 | 0.38 |

Data-only Pareto set (this case and declared settings only): svd_baranger_uncentered, svd_long_uncentered, svd_baranger_centered, svd_long_centered, axis2_centered

## Downstream Doppler and directional WSS

| Method | Velocity NRMSE | Near-wall NRMSE | WSS NRMSE |
|---|---:|---:|---:|
| svd_fixed_uncentered | 0.0158 | 0.0319 | 1.35 |
| svd_energy_uncentered | 0.0158 | 0.0319 | 1.35 |
| svd_baranger_uncentered | 0.0171 | 0.0359 | 0.465 |
| svd_long_uncentered | 0.0154 | 0.0302 | 1.21 |
| svd_local_uncentered | 0.0179 | 0.0337 | 0.959 |
| song_local_reconstruction_uncentered | 0.0566 | 0.0912 | 0.433 |
| svd_temporal_frequency_uncentered | 0.0158 | 0.0319 | 1.35 |
| evd_fixed_uncentered | 0.0158 | 0.0319 | 1.35 |
| song_rsvd_uncentered | 0.0158 | 0.0319 | 1.35 |
| song_rsvd_downsampled_uncentered | 0.0158 | 0.0319 | 1.19 |
| huang_evd_median_uncentered | 4.35 | 6.83 | 1.78 |
| huang_evd_ncc_uncentered | 0.0385 | 0.0666 | 1.58 |
| huang_evd_sigmoid_uncentered | 0.0345 | 0.0477 | 1.63 |
| huang_gaussian_iq_extension_uncentered | 0.0363 | 0.0505 | 1.64 |
| aer_temporal_alpha0.5_uncentered | 0.0337 | 0.0516 | 0.146 |
| aer_temporal_alpha0.25_uncentered | 0.377 | 0.42 | 0.803 |
| aer_spatial_alpha0.5_uncentered | 0.0337 | 0.0516 | 0.146 |
| aer_spatial_alpha0.25_uncentered | 0.377 | 0.42 | 0.803 |
| aer_noise_floor_uncentered | 0.266 | 0.346 | 1.59 |
| axis1_uncentered | 0.0178 | 0.0333 | 0.94 |
| axis1_rigid_uncentered | 0.0178 | 0.0333 | 0.94 |
| axis2_uncentered | 0.602 | 0.663 | 1.05 |
| axis2_with_memory_uncentered | 0.628 | 0.681 | 1.04 |
| svd_fixed_centered | 0.0162 | 0.0339 | 1.45 |
| svd_energy_centered | 0.0162 | 0.0339 | 1.45 |
| svd_baranger_centered | 0.0179 | 0.0395 | 0.513 |
| svd_long_centered | 0.0155 | 0.0309 | 1.33 |
| svd_local_centered | 0.0182 | 0.0358 | 1.15 |
| song_local_reconstruction_centered | 0.0577 | 0.0979 | 0.452 |
| svd_temporal_frequency_centered | 0.63 | 0.641 | 0.987 |
| evd_fixed_centered | 0.0162 | 0.0339 | 1.45 |
| song_rsvd_centered | 0.0162 | 0.0339 | 1.45 |
| song_rsvd_downsampled_centered | 0.0162 | 0.0339 | 1.32 |
| huang_evd_median_centered | 4.35 | 6.87 | 1.53 |
| huang_evd_ncc_centered | 0.044 | 0.0922 | 1.53 |
| huang_evd_sigmoid_centered | 0.0365 | 0.0549 | 1.67 |
| huang_gaussian_iq_extension_centered | 0.0386 | 0.0556 | 1.69 |
| aer_temporal_alpha0.5_centered | 0.0331 | 0.0509 | 0.146 |
| aer_temporal_alpha0.25_centered | 0.374 | 0.418 | 0.801 |
| aer_spatial_alpha0.5_centered | 0.0331 | 0.0509 | 0.146 |
| aer_spatial_alpha0.25_centered | 0.374 | 0.418 | 0.801 |
| aer_noise_floor_centered | 0.274 | 0.346 | 1.63 |
| axis1_centered | 0.0332 | 0.0877 | 1.65 |
| axis1_rigid_centered | 0.0332 | 0.0877 | 1.65 |
| axis2_centered | 0.026 | 0.0561 | 0.364 |
| axis2_with_memory_centered | 0.0285 | 0.0591 | 0.351 |
| ica_circular_complex | 0.0281 | 0.0611 | 1.66 |
| rpca_pcp | 0.15 | 0.172 | 1.53 |
| butterworth_hp_0.01_prf | 0.706 | 0.664 | 0.971 |
| butterworth_hp_0.03_prf | 0.0443 | 0.109 | 0.3 |
| butterworth_hp_0.08_prf | 0.387 | 1.16 | 2.56 |
| butterworth_hp_causal_0.03_prf | 0.082 | 0.12 | 0.286 |
| oracle_svd_preservation_centeredFalse | 0.0158 | 0.0319 | 1.35 |
| oracle_svd_preservation_centeredTrue | 0.0162 | 0.0339 | 1.45 |
| oracle_aer_iq | 0.0322 | 0.0482 | 0.665 |
