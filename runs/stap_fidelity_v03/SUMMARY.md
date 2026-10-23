# First angle-time STAP study

Frozen statistical-IQ study. All methods retain the same views. Calibration is independently acquired.
Oracle SVD uses same-case blood truth and is evaluation-only. LCMV is a separate single-target detector.

| Method | Completed | Median blood IQ distortion | Median tissue retention | Median spectral TV | Median noise gain |
|---|---:|---:|---:|---:|---:|
| unfiltered | 42/42 | 0 | 1 | 0 | 1 |
| butterworth_30Hz | 42/42 | 0.2904 | 3.287e-05 | 0.06922 | 0.9664 |
| svd_joint_matched | 42/42 | 0.5069 | 3.009e-07 | 0.3315 | 0.8565 |
| svd_per_view_matched | 42/42 | 0.5069 | 2.134e-07 | 0.3315 | 0.8565 |
| svd_joint_long | 42/42 | 0.0867 | 2.364e-05 | 0.007382 | 0.9835 |
| baranger_joint_long | 42/42 | 0.2151 | 2.244e-05 | 0.045 | 0.3225 |
| calibrated_demodulated_svd | 42/42 | 0.5069 | 2.307e-07 | 0.3315 | 0.8565 |
| stap_orthogonal | 42/42 | 0.4852 | 6.448e-06 | 0.298 | 0.9526 |
| stap_covariance | 42/42 | 0.9639 | 0.0006554 | 0.2796 | 1.393 |
| stap_covariance_16taps | 42/42 | 0.871 | 0.002271 | 0.2884 | 1.123 |
| stap_training_starved_control | 42/42 | 0.5294 | 0.000532 | 0.3129 | 1.309 |
| oracle_svd_joint_matched | 42/42 | 0.6007 | 3.05e-07 | 0.3268 | 0.0892 |

Medians pool deliberately different conditions, including failure controls. Use the condition records for interpretation.
Velocity spread, acceleration, calibration error, equal blood/tissue velocity and insufficient training are explicit stress cases.
