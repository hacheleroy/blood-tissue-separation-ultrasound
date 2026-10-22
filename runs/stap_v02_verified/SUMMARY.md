# First angle-time STAP study

Frozen statistical-IQ study. All methods retain the same views. Calibration is independently acquired.
Oracle SVD uses same-case blood truth and is evaluation-only. LCMV is a separate single-target detector.

| Method | Completed | Median blood retention | Median tissue retention | Median axial NRMSE | Median WSS NRMSE |
|---|---:|---:|---:|---:|---:|
| unfiltered | 42/42 | 1 | 1 | 0.7656 | 1.469 |
| butterworth_30Hz | 42/42 | 0.9123 | 3.287e-05 | 0.1784 | 0.358 |
| svd_joint_matched | 42/42 | 0.5063 | 3.009e-07 | 0.09686 | 0.7055 |
| svd_per_view_matched | 42/42 | 0.5063 | 2.134e-07 | 0.0972 | 0.7066 |
| svd_joint_long | 42/42 | 0.9933 | 2.364e-05 | 0.07979 | 0.2635 |
| baranger_joint_long | 42/42 | 0.9426 | 2.244e-05 | 0.02548 | 0.1958 |
| calibrated_demodulated_svd | 42/42 | 0.5063 | 2.307e-07 | 0.09709 | 0.707 |
| stap_orthogonal | 42/42 | 0.5731 | 6.448e-06 | 0.4042 | 1.076 |
| stap_covariance | 42/42 | 0.7745 | 0.0006554 | 0.5716 | 1.184 |
| stap_covariance_16taps | 42/42 | 1.246 | 0.002271 | 0.4496 | 1.155 |
| stap_training_starved_control | 42/42 | 0.577 | 0.000532 | 0.4042 | 1.118 |
| oracle_svd_joint_matched | 42/42 | 0.3671 | 3.05e-07 | 0.2075 | 0.9764 |

Medians pool deliberately different conditions, including failure controls. Use the condition records for interpretation.
Velocity spread, acceleration, calibration error, equal blood/tissue velocity and insufficient training are explicit stress cases.
