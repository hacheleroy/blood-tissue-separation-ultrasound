# Benchmark run

Component metrics use the mixture-fitted operator, with no fitted gain or phase correction.
Experimental anatomy and processed velocity maps do not supply independent component truth.

| Method | Status | Blood IQ distortion | Tissue power retention | Blood power retention | Spectral TV | Noise gain | Time (s) |
|---|---|---:|---:|---:|---:|---:|---:|
| unfiltered | ok | — | — | — | — | — | 0.00502 |
| svd_fixed_uncentered | ok | — | — | — | — | — | 0.0146 |
| svd_energy_uncentered | ok | — | — | — | — | — | 0.014 |
| svd_baranger_uncentered | ok | — | — | — | — | — | 0.0162 |
| svd_long_uncentered | ok | — | — | — | — | — | 0.0162 |
| svd_local_uncentered | ok | — | — | — | — | — | 0.136 |
| song_local_reconstruction_uncentered | ok | — | — | — | — | — | 2.91 |
| svd_temporal_frequency_uncentered | ok | — | — | — | — | — | 0.0137 |
| evd_fixed_uncentered | ok | — | — | — | — | — | 0.0111 |
| song_rsvd_uncentered | ok | — | — | — | — | — | 0.00949 |
| song_rsvd_downsampled_uncentered | ok | — | — | — | — | — | 0.00863 |
| huang_evd_median_uncentered | ok | — | — | — | — | — | 0.0103 |
| huang_evd_ncc_uncentered | ok | — | — | — | — | — | 0.0125 |
| huang_evd_sigmoid_uncentered | ok | — | — | — | — | — | 0.013 |
| huang_gaussian_iq_extension_uncentered | ok | — | — | — | — | — | 0.045 |
| aer_temporal_alpha0.5_uncentered | ok | — | — | — | — | — | 0.00885 |
| aer_temporal_alpha0.25_uncentered | ok | — | — | — | — | — | 0.00861 |
| aer_spatial_alpha0.5_uncentered | ok | — | — | — | — | — | 18 |
| aer_spatial_alpha0.25_uncentered | ok | — | — | — | — | — | 16.7 |
| aer_noise_floor_uncentered | ok | — | — | — | — | — | 0.00773 |
| axis1_uncentered | ok | — | — | — | — | — | 0.169 |
| axis1_rigid_uncentered | ok | — | — | — | — | — | 0.177 |
| axis2_uncentered | ok | — | — | — | — | — | 6.59 |
| axis2_with_memory_uncentered | ok | — | — | — | — | — | 7.84 |
| svd_fixed_centered | ok | — | — | — | — | — | 0.0103 |
| svd_energy_centered | ok | — | — | — | — | — | 0.00987 |
| svd_baranger_centered | ok | — | — | — | — | — | 0.0111 |
| svd_long_centered | ok | — | — | — | — | — | 0.0108 |
| svd_local_centered | ok | — | — | — | — | — | 0.0951 |
| song_local_reconstruction_centered | ok | — | — | — | — | — | 2.09 |
| svd_temporal_frequency_centered | ok | — | — | — | — | — | 0.0115 |
| evd_fixed_centered | ok | — | — | — | — | — | 0.00662 |
| song_rsvd_centered | ok | — | — | — | — | — | 0.00659 |
| song_rsvd_downsampled_centered | ok | — | — | — | — | — | 0.00519 |
| huang_evd_median_centered | ok | — | — | — | — | — | 0.00785 |
| huang_evd_ncc_centered | ok | — | — | — | — | — | 0.00824 |
| huang_evd_sigmoid_centered | ok | — | — | — | — | — | 0.00892 |
| huang_gaussian_iq_extension_centered | ok | — | — | — | — | — | 0.0291 |
| aer_temporal_alpha0.5_centered | ok | — | — | — | — | — | 0.00649 |
| aer_temporal_alpha0.25_centered | ok | — | — | — | — | — | 0.00684 |
| aer_spatial_alpha0.5_centered | ok | — | — | — | — | — | 16.6 |
| aer_spatial_alpha0.25_centered | ok | — | — | — | — | — | 15.9 |
| aer_noise_floor_centered | ok | — | — | — | — | — | 0.00807 |
| axis1_centered | ok | — | — | — | — | — | 0.136 |
| axis1_rigid_centered | ok | — | — | — | — | — | 0.228 |
| axis2_centered | ok | — | — | — | — | — | 5.89 |
| axis2_with_memory_centered | ok | — | — | — | — | — | 6.88 |
| ica_circular_complex | ok | — | — | — | — | — | 9.72 |
| rpca_pcp | ok | — | — | — | — | — | 0.374 |
| butterworth_hp_0.01_prf | ok | — | — | — | — | — | 0.0106 |
| butterworth_hp_0.03_prf | ok | — | — | — | — | — | 0.0161 |
| butterworth_hp_0.08_prf | ok | — | — | — | — | — | 0.015 |
| butterworth_hp_causal_0.03_prf | ok | — | — | — | — | — | 0.00913 |
| view_svd_matched | ok | — | — | — | — | — | 0.0338 |
| view_svd_long | ok | — | — | — | — | — | 0.0347 |
| calibrated_view_svd_matched | ok | — | — | — | — | — | 0.0281 |
| stap_orthogonal | ok | — | — | — | — | — | 0.146 |
| stap_covariance | ok | — | — | — | — | — | 5.76 |

Data-only Pareto set (this case and declared settings only):

## Downstream Doppler and directional WSS

| Method | Velocity NRMSE | Near-wall NRMSE | WSS NRMSE |
|---|---:|---:|---:|
