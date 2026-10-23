# Separation fidelity study

Primary comparison: isolated blood complex-IQ distortion versus tissue power retention.
No fitted gain/phase correction. Filter parameters and case bank are frozen before execution.
The equal-static case is an identifiability control; complete separation is not expected.

| Family | Method | Completed | Median blood distortion | Median tissue retention | Median spectral TV | ICA convergence |
|---|---|---:|---:|---:|---:|---:|
| statistical | unfiltered | 15/15 | 0 | 1 | 0 | n/a |
| statistical | svd_fixed_uncentered | 15/15 | 0.2913 | 7.879e-07 | 0.0648 | n/a |
| statistical | svd_fixed_centered | 15/15 | 0.5475 | 7.495e-07 | 0.2723 | n/a |
| statistical | svd_baranger_uncentered | 15/15 | 0.3811 | 5.664e-07 | 0.1187 | n/a |
| statistical | svd_baranger_centered | 15/15 | 0.574 | 5.537e-07 | 0.2719 | n/a |
| statistical | svd_local_uncentered | 15/15 | 0.2994 | 1.174e-05 | 0.0712 | n/a |
| statistical | svd_local_centered | 15/15 | 0.5531 | 1.195e-05 | 0.2716 | n/a |
| statistical | svd_long_uncentered | 15/15 | 0.1388 | 8.259e-07 | 0.01486 | n/a |
| statistical | svd_long_centered | 15/15 | 0.1786 | 8.152e-07 | 0.03446 | n/a |
| statistical | butterworth_hp_0.01_prf | 15/15 | 0.4613 | 0.3266 | 0.1303 | n/a |
| statistical | butterworth_hp_0.03_prf | 15/15 | 0.6357 | 0.000493 | 0.3348 | n/a |
| statistical | butterworth_hp_0.08_prf | 15/15 | 0.9035 | 3.793e-06 | 0.6993 | n/a |
| statistical | aer_temporal_alpha0.5_uncentered | 15/15 | 0.874 | 1.978e-05 | 0.04613 | n/a |
| statistical | aer_temporal_alpha0.25_uncentered | 15/15 | 0.6469 | 0.004446 | 0.03297 | n/a |
| statistical | aer_temporal_alpha0.5_centered | 15/15 | 0.876 | 1.978e-05 | 0.08674 | n/a |
| statistical | aer_temporal_alpha0.25_centered | 15/15 | 0.6555 | 0.004316 | 0.06322 | n/a |
| statistical | aer_spatial_alpha0.5_uncentered | 15/15 | 0.8777 | 1.84e-05 | 0.02344 | n/a |
| statistical | aer_spatial_alpha0.5_centered | 15/15 | 0.8797 | 1.839e-05 | 0.04801 | n/a |
| statistical | axis1_uncentered | 15/15 | 0.1287 | 2.077e-05 | 0.01471 | n/a |
| statistical | axis1_centered | 15/15 | 0.639 | 1.152e-07 | 0.2219 | n/a |
| statistical | axis1_rigid_uncentered | 15/15 | 0.1287 | 2.077e-05 | 0.01482 | n/a |
| statistical | axis1_rigid_centered | 15/15 | 0.6644 | 1.152e-07 | 0.2219 | n/a |
| statistical | axis2_uncentered | 15/15 | 3.361 | 0.2465 | 0.2421 | n/a |
| statistical | axis2_centered | 15/15 | 0.7876 | 6.88e-10 | 0.1123 | n/a |
| statistical | axis2_with_memory_uncentered | 15/15 | 2.679 | 0.2382 | 0.1829 | n/a |
| statistical | axis2_with_memory_centered | 15/15 | 0.758 | 8.513e-10 | 0.1491 | n/a |
| statistical | ica_circular_complex | 15/15 | 0.7238 | 1.112e-07 | 0.3092 | 5/352 |
| statistical | rpca_pcp | 15/15 | n/a | n/a | n/a | n/a |
| statistical | ica_uncentered_retained_300 | 15/15 | 3.286 | 0.5539 | 0.1152 | 43/352 |
| statistical | ica_centered_retained_300 | 15/15 | 0.8565 | 8.323e-10 | 0.0855 | 5/352 |
| statistical | ica_centered_retained_1500 | 15/15 | 0.8564 | 8.362e-10 | 0.1268 | 7/352 |
| statistical | axis1_no_memory | 15/15 | 0.1313 | 4.085e-05 | 0.0143 | n/a |
| statistical | axis1_no_phase | 15/15 | 0.1287 | 2.077e-05 | 0.01471 | n/a |
| statistical | axis2_no_phase | 15/15 | 1.639 | 0.4654 | 0.164 | n/a |
| statistical | axis1_phase_error | 15/15 | 0.1287 | 2.077e-05 | 0.01471 | n/a |
| statistical | axis2_phase_error | 15/15 | 2.236 | 0.2607 | 0.1263 | n/a |
| statistical | oracle_svd_iq | 15/15 | 0.7348 | 6.222e-07 | 0.1459 | n/a |
| statistical | oracle_svd_distortion_025 | 15/15 | 1.317e-15 | 1 | 2.822e-16 | n/a |
| statistical | oracle_svd_distortion_050 | 15/15 | 0.2097 | 1.726e-06 | 0.02863 | n/a |
| psf_transport | unfiltered | 15/15 | 0 | 1 | 0 | n/a |
| psf_transport | svd_fixed_uncentered | 15/15 | 0.6734 | 0.02254 | 0.3614 | n/a |
| psf_transport | svd_fixed_centered | 15/15 | 0.8076 | 0.001356 | 0.4327 | n/a |
| psf_transport | svd_baranger_uncentered | 15/15 | 0.8806 | 1.643e-06 | 0.6071 | n/a |
| psf_transport | svd_baranger_centered | 15/15 | 0.9999 | 5.417e-08 | 0.8832 | n/a |
| psf_transport | svd_local_uncentered | 15/15 | 0.68 | 0.01428 | 0.3678 | n/a |
| psf_transport | svd_local_centered | 15/15 | 0.8276 | 0.0006925 | 0.4459 | n/a |
| psf_transport | svd_long_uncentered | 15/15 | 0.4092 | 0.1434 | 0.1123 | n/a |
| psf_transport | svd_long_centered | 15/15 | 0.4683 | 0.06129 | 0.1208 | n/a |
| psf_transport | butterworth_hp_0.01_prf | 15/15 | 0.5357 | 0.05648 | 0.1299 | n/a |
| psf_transport | butterworth_hp_0.03_prf | 15/15 | 0.7257 | 0.0001946 | 0.4026 | n/a |
| psf_transport | butterworth_hp_0.08_prf | 15/15 | 0.9901 | 1.125e-06 | 0.7031 | n/a |
| psf_transport | aer_temporal_alpha0.5_uncentered | 15/15 | 0.9161 | 8.287e-05 | 0.03151 | n/a |
| psf_transport | aer_temporal_alpha0.25_uncentered | 15/15 | 0.7465 | 0.006993 | 0.02184 | n/a |
| psf_transport | aer_temporal_alpha0.5_centered | 15/15 | 0.9194 | 7.265e-05 | 0.346 | n/a |
| psf_transport | aer_temporal_alpha0.25_centered | 15/15 | 0.7612 | 0.005532 | 0.1033 | n/a |
| psf_transport | aer_spatial_alpha0.5_uncentered | 15/15 | 0.8575 | 8.254e-05 | 0.02267 | n/a |
| psf_transport | aer_spatial_alpha0.5_centered | 15/15 | 0.8633 | 7.462e-05 | 0.07997 | n/a |
| psf_transport | axis1_uncentered | 15/15 | 0.1036 | 0.04415 | 0.005689 | n/a |
| psf_transport | axis1_centered | 15/15 | 0.5959 | 0.006263 | 0.1495 | n/a |
| psf_transport | axis1_rigid_uncentered | 15/15 | 0.1036 | 0.04415 | 0.005689 | n/a |
| psf_transport | axis1_rigid_centered | 15/15 | 0.5959 | 0.006263 | 0.1495 | n/a |
| psf_transport | axis2_uncentered | 15/15 | 5.2 | 0.2885 | 0.2174 | n/a |
| psf_transport | axis2_centered | 15/15 | 1.087 | 0.005083 | 0.4185 | n/a |
| psf_transport | axis2_with_memory_uncentered | 15/15 | 3.642 | 0.3086 | 0.2155 | n/a |
| psf_transport | axis2_with_memory_centered | 15/15 | 1.003 | 0.03603 | 0.3065 | n/a |
| psf_transport | ica_circular_complex | 15/15 | 0.7977 | 0.0005028 | 0.3768 | 14/360 |
| psf_transport | rpca_pcp | 15/15 | n/a | n/a | n/a | n/a |
| psf_transport | ica_uncentered_retained_300 | 15/15 | 9.298 | 0.4435 | 0.2715 | 20/360 |
| psf_transport | ica_centered_retained_300 | 15/15 | 1.162 | 0.006284 | 0.4186 | 14/360 |
| psf_transport | ica_centered_retained_1500 | 15/15 | 1.222 | 0.00613 | 0.4493 | 15/360 |
| psf_transport | axis1_no_memory | 15/15 | 0.1146 | 0.0151 | 0.005585 | n/a |
| psf_transport | axis1_no_phase | 15/15 | 0.1036 | 0.04415 | 0.005689 | n/a |
| psf_transport | axis2_no_phase | 15/15 | 3.389 | 0.4633 | 0.2123 | n/a |
| psf_transport | axis1_phase_error | 15/15 | 0.1036 | 0.04415 | 0.005689 | n/a |
| psf_transport | axis2_phase_error | 15/15 | 3.653 | 0.3339 | 0.213 | n/a |
| psf_transport | oracle_svd_iq | 15/15 | 1 | 0 | 0.161 | n/a |
| psf_transport | oracle_svd_distortion_025 | 15/15 | 8.942e-16 | 1 | 1.657e-16 | n/a |
| psf_transport | oracle_svd_distortion_050 | 15/15 | 8.942e-16 | 1 | 2.359e-16 | n/a |

Use per-condition Pareto sets; pooled medians are descriptive, not a universal ranking.
Infeasible oracle searches retain their full candidate curves. RPCA uses labeled refitted injection responses.
The 2-D PSF model includes transport and boundary mixing but omits channel propagation and RBC physics.
