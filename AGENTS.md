# Project working agreement

The user authorized autonomous validation and direct pushes to `main` on
23 October 2026. Work on `main`; do not create a development branch or pull
request, or ask for routine external validation, unless the user requests it.

Before pushing, run numerical checks appropriate to the change, verify frozen
study archives and source identity, and inspect the final diff. Publish with an
expected-head lease and no force; reconcile concurrent changes before retrying.
Check the resulting GitHub Actions run and fix failures autonomously.

Preserve historical results/configurations and their hashes. Develop or tune on
declared development/training data. Freeze numerical code, configuration and
selection rules before inspecting a new held-out bank; never relabel seen
evaluation data as held out. Keep known motion and source truth out of regular
filtering and gating APIs.

Retain the research-note names:

- Axis 1: Motion-compensated virtual long ensemble
- Axis 2: Motion-demodulated multi-lag separation
- Axis 3: Local ultrasound STAP
