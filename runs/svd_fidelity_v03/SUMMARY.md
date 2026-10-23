# Literature SVD/EVD comparison

Frozen statistical-IQ cases. See SVD_LITERATURE_AUDIT.md for source fidelity and provisional arms.

| Method | Completed | Median blood IQ distortion | Median tissue retention | Median blood retention | Median spectral TV |
|---|---:|---:|---:|---:|---:|
| svd_fixed_uncentered | 24/24 | 0.4588 | 2.26e-06 | 0.6797 | 0.2709 |
| svd_fixed_centered | 24/24 | 0.5573 | 1.563e-06 | 0.5592 | 0.3879 |
| svd_energy_uncentered | 24/24 | 0.4588 | 2.26e-06 | 0.6797 | 0.2709 |
| svd_baranger_uncentered | 24/24 | 0.4701 | 1.768e-06 | 0.6745 | 0.2701 |
| svd_baranger_centered | 24/24 | 0.5641 | 1.144e-06 | 0.5548 | 0.3871 |
| svd_local_uncentered | 24/24 | 0.4585 | 2.787e-05 | 0.6796 | 0.2714 |
| svd_long_uncentered | 24/24 | 0.2369 | 2.423e-06 | 0.9111 | 0.07073 |
| song_local_reconstruction_uncentered | 24/24 | 0.8638 | 8.715e-07 | 0.1123 | 0.5485 |
| svd_temporal_frequency_uncentered | 24/24 | 0.3122 | 0.3333 | 0.8531 | 0.04462 |
| evd_fixed_uncentered | 24/24 | 0.4588 | 2.26e-06 | 0.6797 | 0.2709 |
| song_rsvd_uncentered | 24/24 | 0.4588 | 2.26e-06 | 0.6797 | 0.2709 |
| song_rsvd_downsampled_uncentered | 24/24 | 0.4588 | 3.24e-06 | 0.6797 | 0.2709 |
| huang_evd_median_uncentered | 24/24 | 0.9965 | 4.566e-07 | 0.007084 | 0.9518 |
| huang_evd_ncc_uncentered | 24/24 | 0.9921 | 5.749e-07 | 0.01569 | 0.9329 |
| huang_evd_sigmoid_uncentered | 24/24 | 0.9884 | 5.936e-07 | 0.01634 | 0.9236 |
| huang_gaussian_iq_extension_uncentered | 24/24 | 0.988 | 5.959e-07 | 0.01453 | 0.9071 |
| butterworth_hp_0.01_prf | 24/24 | 0.5007 | 0.3043 | 0.8802 | 0.1405 |
| butterworth_hp_0.03_prf | 24/24 | 0.6031 | 0.004112 | 0.5561 | 0.4044 |
| butterworth_hp_0.08_prf | 24/24 | 0.714 | 1.353e-05 | 0.3264 | 0.6037 |
| butterworth_hp_causal_0.03_prf | 24/24 | 0.893 | 0.01535 | 0.5529 | 0.2985 |

Every failed or skipped method is retained in results.json. Pooled medians do not establish a universal ranking.
Gaussian-localized IQ overlap-add is an extension; paper power aggregation is reported separately.
Song thresholds are a documented reconstruction, pending original-paper/software parity.
