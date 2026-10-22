# Literature SVD/EVD comparison

Frozen statistical-IQ cases. See SVD_LITERATURE_AUDIT.md for source fidelity and provisional arms.

| Method | Completed | Median velocity NRMSE | Median WSS NRMSE | Median blood retention | Median tissue retention |
|---|---:|---:|---:|---:|---:|
| svd_fixed_uncentered | 24/24 | 0.2504 | 21.77 | 0.6797 | 2.26e-06 |
| svd_fixed_centered | 24/24 | 0.4705 | 24.03 | 0.5592 | 1.563e-06 |
| svd_energy_uncentered | 24/24 | 0.2504 | 21.77 | 0.6797 | 2.26e-06 |
| svd_baranger_uncentered | 24/24 | 0.2399 | 3.574 | 0.6745 | 1.768e-06 |
| svd_baranger_centered | 24/24 | 0.4699 | 3.86 | 0.5548 | 1.144e-06 |
| svd_local_uncentered | 24/24 | 0.2382 | 22.52 | 0.6796 | 2.787e-05 |
| svd_long_uncentered | 24/24 | 0.1394 | 20.14 | 0.9111 | 2.423e-06 |
| song_local_reconstruction_uncentered | 24/24 | 0.6294 | 3.443 | 0.1123 | 8.715e-07 |
| svd_temporal_frequency_uncentered | 24/24 | 0.4397 | 3.939 | 0.8531 | 0.3333 |
| evd_fixed_uncentered | 24/24 | 0.2504 | 21.77 | 0.6797 | 2.26e-06 |
| song_rsvd_uncentered | 24/24 | 0.2504 | 21.77 | 0.6797 | 2.26e-06 |
| song_rsvd_downsampled_uncentered | 24/24 | 0.2505 | 22.36 | 0.6797 | 3.24e-06 |
| huang_evd_median_uncentered | 24/24 | 17.03 | 19.01 | 0.007084 | 4.566e-07 |
| huang_evd_ncc_uncentered | 24/24 | 7.932 | 32.36 | 0.01569 | 5.749e-07 |
| huang_evd_sigmoid_uncentered | 24/24 | 7.043 | 32.92 | 0.01634 | 5.936e-07 |
| huang_gaussian_iq_extension_uncentered | 24/24 | 11.19 | 18.97 | 0.01453 | 5.959e-07 |
| butterworth_hp_0.01_prf | 24/24 | 0.3121 | 0.9959 | 0.8802 | 0.3043 |
| butterworth_hp_0.03_prf | 24/24 | 0.5745 | 1.655 | 0.5561 | 0.004112 |
| butterworth_hp_0.08_prf | 24/24 | 2.484 | 12.11 | 0.3264 | 1.353e-05 |
| butterworth_hp_causal_0.03_prf | 24/24 | 0.4525 | 1.077 | 0.5529 | 0.01535 |

Every failed or skipped method is retained in results.json. Pooled medians do not establish a universal ranking.
Gaussian-localized IQ overlap-add is an extension; paper power aggregation is reported separately.
Song thresholds are a documented reconstruction, pending original-paper/software parity.
