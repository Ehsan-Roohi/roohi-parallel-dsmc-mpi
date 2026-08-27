# Parallel DSMC Solver with MPI

A legacy MPI-parallel Direct Simulation Monte Carlo solver associated with the accepted parallel-DSMC journal project.

## Pre-release status

This repository was reviewed and approved for public release by Ehsan Roohi on 2026-08-27. The original working files remain preserved separately.

## Contents

- `src/`: 1 conservatively selected research-code files.
- `docs/SOURCE_MANIFEST.csv`: exact mapping to the original local files, including SHA-256 hashes.
- `docs/EXCLUDED_FILES.csv`: duplicate, attributed-to-others, private, malformed, or otherwise withheld files.
- `docs/DATA_IO_REFERENCES.csv`: locations that read or write data and therefore need input/output documentation.
- `OWNERSHIP_REVIEW.md`: release gate and authorship checklist.

## Languages

C, MPI

## License

Released under the MIT License. See `LICENSE`.

## Notes

The located parallel implementation is C/MPI rather than Fortran. No separate MPI/OpenMP Fortran source was found in the workspace. Compilation requires an MPI C compiler such as `mpicc`; scientific verification and input documentation remain required.
