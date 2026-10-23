# Benchmark run

Component metrics use the mixture-fitted operator, with no fitted gain or phase correction.
Experimental anatomy and processed velocity maps do not supply independent component truth.

| Method | Status | Blood IQ distortion | Tissue power retention | Blood power retention | Spectral TV | Noise gain | Time (s) |
|---|---|---:|---:|---:|---:|---:|---:|
| unfiltered | ok | — | — | — | — | — | 0.0395 |
| svd_w32_k1 | ok | — | — | — | — | — | 0.118 |
| svd_w32_k4 | ok | — | — | — | — | — | 0.109 |
| svd_w256_k1 | ok | — | — | — | — | — | 0.367 |
| svd_w256_k4 | ok | — | — | — | — | — | 0.365 |
| svd_w768_k1 | ok | — | — | — | — | — | 0.961 |
| svd_w768_k4 | ok | — | — | — | — | — | 0.943 |
| baranger_w256 | ok | — | — | — | — | — | 0.387 |
| maresca_k26_w60 | ok | — | — | — | — | — | 0.198 |
| song_pdf_paper | ok | — | — | — | — | — | 0.398 |
| song_pdf_none | ok | — | — | — | — | — | 0.399 |
| lee_temporal_w256 | ok | — | — | — | — | — | 0.224 |
| lee_temporal_w768 | ok | — | — | — | — | — | 0.778 |
| lee_temporal_absolute26 | ok | — | — | — | — | — | 0.751 |
| lee_temporal_unit_extension | ok | — | — | — | — | — | 0.751 |
| lee_spatial_w256 | ok | — | — | — | — | — | 1.69 |
| lee_spatial_w768 | ok | — | — | — | — | — | 0.877 |
| lee_spatial_absolute26 | ok | — | — | — | — | — | 0.86 |
| lee_spatial_unit_extension | ok | — | — | — | — | — | 0.882 |
| axis1_memory0_rigid0 | ok | — | — | — | — | — | 0.689 |
| axis1_memory1_rigid0 | ok | — | — | — | — | — | 0.572 |
| axis1_memory1_rigid1 | ok | — | — | — | — | — | 1.23 |
| axis2_c0_guard0 | ok | — | — | — | — | — | 14.9 |
| axis2_c0_guard1 | ok | — | — | — | — | — | 12.9 |
| axis2_c1_guard0 | ok | — | — | — | — | — | 16.1 |
| axis2_c1_guard1 | ok | — | — | — | — | — | 9.66 |
| view_svd_matched | ok | — | — | — | — | — | 0.239 |
| stap_orthogonal | ok | — | — | — | — | — | 0.745 |
| stap_covariance | ok | — | — | — | — | — | 36.7 |

Data-only Pareto set (this case and declared settings only):

## Downstream Doppler and directional WSS

| Method | Velocity NRMSE | Near-wall NRMSE | WSS NRMSE |
|---|---:|---:|---:|
