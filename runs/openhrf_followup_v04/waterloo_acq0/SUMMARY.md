# Benchmark run

Component metrics use the mixture-fitted operator, with no fitted gain or phase correction.
Experimental anatomy and processed velocity maps do not supply independent component truth.

| Method | Status | Blood IQ distortion | Tissue power retention | Blood power retention | Spectral TV | Noise gain | Time (s) |
|---|---|---:|---:|---:|---:|---:|---:|
| unfiltered | ok | — | — | — | — | — | 0.0366 |
| svd_w32_k1 | ok | — | — | — | — | — | 0.081 |
| svd_w32_k4 | ok | — | — | — | — | — | 0.0872 |
| svd_w256_k1 | ok | — | — | — | — | — | 0.271 |
| svd_w256_k4 | ok | — | — | — | — | — | 0.254 |
| svd_w768_k1 | ok | — | — | — | — | — | 0.748 |
| svd_w768_k4 | ok | — | — | — | — | — | 0.778 |
| baranger_w256 | ok | — | — | — | — | — | 0.275 |
| maresca_k26_w60 | ok | — | — | — | — | — | 0.127 |
| song_pdf_paper | ok | — | — | — | — | — | 0.265 |
| song_pdf_none | ok | — | — | — | — | — | 0.263 |
| lee_temporal_w256 | ok | — | — | — | — | — | 0.15 |
| lee_temporal_w768 | ok | — | — | — | — | — | 0.516 |
| lee_temporal_absolute26 | ok | — | — | — | — | — | 0.507 |
| lee_temporal_unit_extension | ok | — | — | — | — | — | 0.692 |
| lee_spatial_w256 | ok | — | — | — | — | — | 1.51 |
| lee_spatial_w768 | ok | — | — | — | — | — | 0.856 |
| lee_spatial_absolute26 | ok | — | — | — | — | — | 0.874 |
| lee_spatial_unit_extension | ok | — | — | — | — | — | 0.826 |
| axis1_memory0_rigid0 | ok | — | — | — | — | — | 0.588 |
| axis1_memory1_rigid0 | ok | — | — | — | — | — | 0.597 |
| axis1_memory1_rigid1 | ok | — | — | — | — | — | 1.12 |
| axis2_c0_guard0 | ok | — | — | — | — | — | 10.5 |
| axis2_c0_guard1 | ok | — | — | — | — | — | 5.73 |
| axis2_c1_guard0 | ok | — | — | — | — | — | 8.14 |
| axis2_c1_guard1 | ok | — | — | — | — | — | 8.06 |
| view_svd_matched | ok | — | — | — | — | — | 0.218 |
| stap_orthogonal | ok | — | — | — | — | — | 0.664 |
| stap_covariance | ok | — | — | — | — | — | 28.5 |

Data-only Pareto set (this case and declared settings only):

## Downstream Doppler and directional WSS

| Method | Velocity NRMSE | Near-wall NRMSE | WSS NRMSE |
|---|---:|---:|---:|
