# Public filenames and internal revision provenance

The GitHub release uses stable, descriptive filenames. During manuscript
revision the same scripts carried internal version numbers to protect completed
runs and preserve checkpoints. The public mapping is:

| Public filename | Internal revision source |
|---|---|
| `PACT_Benchmark.m` | `PACT_Benchmark_v93_reviewer_revision.m` |
| `PACT_SoundSpeedMismatch.m` | `PACT_Benchmark_v94_mismatch_only.m` |
| `PACT_NoiseRobustness.m` | `PACT_Benchmark_v95_noise_only.m` |
| `PACT_AdaptiveCohort.m` | `PACT_Benchmark_v96_adaptive_11phantom_only.m` |
| `PACT_PointTargetResolution.m` | `PACT_Benchmark_v97_resolution_9point_resume.m` |

Analysis files were renamed similarly to stable descriptive names. Numerical
settings and algorithms were not changed by this public renaming.
