# PACT 2-D Direct-Beamforming Benchmark

MATLAB code supporting the revised manuscript:

**A Physics-Based Two-Dimensional Numerical Benchmark of Direct Beamforming
Methods for Photoacoustic Computed Tomography in Heterogeneous Media**

## Repository layout

- `benchmark/PACT_Benchmark.m` — main benchmark, including phantom generation,
  heterogeneous optical/acoustic model construction, finite-aperture detector
  model, k-Wave forward simulation, beamformers, metrics, and exports.
- `benchmark/PACT_SoundSpeedMismatch.m` — sound-speed/delay robustness.
- `benchmark/PACT_NoiseRobustness.m` — alternative-noise stress test.
- `benchmark/PACT_AdaptiveCohort.m` — limited-view adaptive-cohort completion.
- `benchmark/PACT_PointTargetResolution.m` — nine-position point-target study.
- `analysis/` — manuscript post-processing and verification routines.
- `examples/` — small launcher scripts.
- `config/` — archived configuration and run provenance.
- `environment/` — MATLAB/k-Wave environment capture.

## Naming

Public filenames are intentionally stable and descriptive; manuscript-development
version numbers are documented in `VERSION_HISTORY.md` rather than embedded in
the MATLAB filenames.

## Requirements

MATLAB with k-Wave and the MATLAB functions/toolboxes checked by
`PACT_Benchmark.m`.

Before archiving a release, run:

```matlab
environment/PACT_collect_environment_info
```

from the repository root and include the generated environment report.

## Reproducibility data

Processed CSV data supporting the principal manuscript figures and tables should
be deposited in a separate `processed_data/` directory in the final repository
release.
