# Phase-6 development-only spatial registration probe

After observing that the declared local-phase grid could not meet the rigid
leakage bound on **design** cloud 2309, this additional model was introduced
within the design split. Seeds 2333, 2401 and 2411 are still unopened.
The [new context](../examples/esdiva_phase6_unit_shift.py) estimates an axial
and lateral translation separately for each observed frame and view from
the 32×16 tissue ROI. The first-half calibration fixes a tissue magnitude
reference; 127 checkerboard tissue pixels fit a bounded translation and 65
disjoint pixels accept it only if their normalized residual improves by
more than 1%. Search bounds are ±1.25 axial and ±0.75 lateral pixels.
Neither isolated source nor prescribed displacement is passed to the fit.

The accepted translation is applied on the original eight-row IQ tile by a
cyclic Fourier shift and exactly reversed after filtering. This makes the
algebraic roundtrip numerically reversible but **assumes periodic boundaries
on an eight-row tile**. A good roundtrip alone is insufficient: wraparound
may change blood response or introduce artificial tissue agreement. The
[development scorer](../examples/esdiva_phase6_score.py) therefore retains
the source-resolved paired blood NRMSE, separately refitted tissue leakage,
wide/core input audit and no-motion control. Shift only and shift plus the
fixed local-phase field are explicit *exploratory* variants. They cannot be
chosen from the old held-out seeds or reported as a validated physical
registration method.

The original plan's eight phase fraction/width settings remain the bounded
comparison. A spatial candidate can advance to seed 2333 only if **both**
design clouds satisfy the original all-scene paired thresholds and exact
roundtrip, with its settings and code frozen first. If it fails, preserve the
negative result and leave the new reserve unopened. The physical bridge does
not inherit a cyclic core boundary assumption without a separate preacquisition
port and equivalence/edge validation.
