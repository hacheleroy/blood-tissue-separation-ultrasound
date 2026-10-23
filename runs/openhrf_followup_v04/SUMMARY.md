# Longer Waterloo experimental processing pilot

Two distinct subjects; 1024 reconstructed frames each; first 256 calibrate nuisance ROI/velocity, last 768 evaluate.
All parameters fixed before evaluation. No independent source/velocity/phase truth, no oracle and no verified contrast status.
Display panels are normalized separately; their appearance is not a separation-accuracy score.

| Acquisition / subject | Method | Status | Bright-ROI output/input power | Median lag-one coherence | Maximum local operator gain |
|---|---|---|---:|---:|---:|
| waterloo_acq0 / 1 | unfiltered | ok | 1 | 0.999826 | n/a |
| waterloo_acq0 / 1 | svd_w32_k1 | ok | 0.00269006 | 0.841643 | n/a |
| waterloo_acq0 / 1 | svd_w32_k4 | ok | 4.6601e-06 | 0.174218 | n/a |
| waterloo_acq0 / 1 | svd_w256_k1 | ok | 0.065035 | 0.995782 | n/a |
| waterloo_acq0 / 1 | svd_w256_k4 | ok | 0.00113506 | 0.964478 | n/a |
| waterloo_acq0 / 1 | svd_w768_k1 | ok | 0.107202 | 0.999597 | n/a |
| waterloo_acq0 / 1 | svd_w768_k4 | ok | 0.013035 | 0.998166 | n/a |
| waterloo_acq0 / 1 | baranger_w256 | ok | 3.98871e-06 | 0.593275 | n/a |
| waterloo_acq0 / 1 | maresca_k26_w60 | ok | 7.03494e-07 | 0.522851 | n/a |
| waterloo_acq0 / 1 | song_pdf_paper | ok | 7.78831e-13 | 0.414569 | n/a |
| waterloo_acq0 / 1 | song_pdf_none | ok | 7.21662e-06 | 0.446813 | n/a |
| waterloo_acq0 / 1 | lee_temporal_w256 | ok | 6.79834e-10 | 0.924381 | n/a |
| waterloo_acq0 / 1 | lee_temporal_w768 | ok | 5.0864e-10 | 0.975706 | n/a |
| waterloo_acq0 / 1 | lee_temporal_absolute26 | ok | 7.61463e-09 | 0.229958 | n/a |
| waterloo_acq0 / 1 | lee_temporal_unit_extension | ok | 0.000800524 | 0.975706 | n/a |
| waterloo_acq0 / 1 | lee_spatial_w256 | ok | 2.25726e-10 | 0.924381 | n/a |
| waterloo_acq0 / 1 | lee_spatial_w768 | ok | 5.07977e-10 | 0.975706 | n/a |
| waterloo_acq0 / 1 | lee_spatial_absolute26 | ok | 7.60646e-09 | 0.229979 | n/a |
| waterloo_acq0 / 1 | lee_spatial_unit_extension | ok | 0.000800524 | 0.975706 | n/a |
| waterloo_acq0 / 1 | axis1_memory0_rigid0 | ok | 0.000621117 | 0.837154 | 1 |
| waterloo_acq0 / 1 | axis1_memory1_rigid0 | ok | 0.0165728 | 0.998121 | 1 |
| waterloo_acq0 / 1 | axis1_memory1_rigid1 | ok | 0.0165728 | 0.998121 | 1 |
| waterloo_acq0 / 1 | axis2_c0_guard0 | ok | 0.216269 | 0.999389 | 55.3917 |
| waterloo_acq0 / 1 | axis2_c0_guard1 | ok | 0.195241 | 0.999712 | 1 |
| waterloo_acq0 / 1 | axis2_c1_guard0 | ok | 0.0240305 | 0.998124 | 10.7741 |
| waterloo_acq0 / 1 | axis2_c1_guard1 | ok | 0.0222184 | 0.998573 | 1 |
| waterloo_acq0 / 1 | view_svd_matched | ok | 3.04099e-06 | 0.413587 | n/a |
| waterloo_acq0 / 1 | stap_orthogonal | ok | 0.000301451 | 0.803837 | n/a |
| waterloo_acq0 / 1 | stap_covariance | ok | 0.0753069 | 0.993562 | n/a |
| waterloo_acq12 / 3 | unfiltered | ok | 1 | 0.999889 | n/a |
| waterloo_acq12 / 3 | svd_w32_k1 | ok | 0.00337425 | 0.851572 | n/a |
| waterloo_acq12 / 3 | svd_w32_k4 | ok | 6.90668e-06 | 0.241021 | n/a |
| waterloo_acq12 / 3 | svd_w256_k1 | ok | 0.0717501 | 0.995686 | n/a |
| waterloo_acq12 / 3 | svd_w256_k4 | ok | 0.00193804 | 0.969843 | n/a |
| waterloo_acq12 / 3 | svd_w768_k1 | ok | 0.13119 | 0.99956 | n/a |
| waterloo_acq12 / 3 | svd_w768_k4 | ok | 0.0177437 | 0.997122 | n/a |
| waterloo_acq12 / 3 | baranger_w256 | ok | 7.40026e-06 | 0.706678 | n/a |
| waterloo_acq12 / 3 | maresca_k26_w60 | ok | 7.88956e-07 | 0.442372 | n/a |
| waterloo_acq12 / 3 | song_pdf_paper | ok | 8.041e-13 | 0.460282 | n/a |
| waterloo_acq12 / 3 | song_pdf_none | ok | 6.34092e-06 | 0.407835 | n/a |
| waterloo_acq12 / 3 | lee_temporal_w256 | ok | 9.94657e-10 | 0.900945 | n/a |
| waterloo_acq12 / 3 | lee_temporal_w768 | ok | 7.18408e-10 | 0.965633 | n/a |
| waterloo_acq12 / 3 | lee_temporal_absolute26 | ok | 9.02372e-09 | 0.25528 | n/a |
| waterloo_acq12 / 3 | lee_temporal_unit_extension | ok | 0.00114642 | 0.965633 | n/a |
| waterloo_acq12 / 3 | lee_spatial_w256 | ok | 3.30257e-10 | 0.900945 | n/a |
| waterloo_acq12 / 3 | lee_spatial_w768 | ok | 7.17473e-10 | 0.965633 | n/a |
| waterloo_acq12 / 3 | lee_spatial_absolute26 | ok | 9.01384e-09 | 0.255302 | n/a |
| waterloo_acq12 / 3 | lee_spatial_unit_extension | ok | 0.00114642 | 0.965633 | n/a |
| waterloo_acq12 / 3 | axis1_memory0_rigid0 | ok | 0.00106518 | 0.850049 | 1 |
| waterloo_acq12 / 3 | axis1_memory1_rigid0 | ok | 0.0235774 | 0.997884 | 1 |
| waterloo_acq12 / 3 | axis1_memory1_rigid1 | ok | 0.0235774 | 0.997884 | 1 |
| waterloo_acq12 / 3 | axis2_c0_guard0 | ok | 0.231179 | 0.999455 | 36.7838 |
| waterloo_acq12 / 3 | axis2_c0_guard1 | ok | 0.215002 | 0.999757 | 1 |
| waterloo_acq12 / 3 | axis2_c1_guard0 | ok | 0.0269869 | 0.997608 | 14.9641 |
| waterloo_acq12 / 3 | axis2_c1_guard1 | ok | 0.0250644 | 0.998088 | 1 |
| waterloo_acq12 / 3 | view_svd_matched | ok | 4.20705e-06 | 0.402602 | n/a |
| waterloo_acq12 / 3 | stap_orthogonal | ok | 0.00833845 | 0.992174 | n/a |
| waterloo_acq12 / 3 | stap_covariance | ok | 0.121343 | 0.998071 | n/a |

The ROI is the brightest calibration-power quartile; its tissue-only status is unverified.
A local cap does not guarantee full adaptive-pipeline stability or correct source identification.
Absolute event timing and both uncompounded views are retained; each evaluation segment lasts about 0.256 s.
