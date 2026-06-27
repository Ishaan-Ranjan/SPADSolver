# SPADSolver

Numerical core for SPAD (Single-Photon Avalanche Diode) solver calculations, with C++ and Python implementations.

## Repository layout

```
SPADSolver/
├── cpp/                 # C++ library and tests
├── python/              # Python package and pytest suite
└── shared/              # Constants and cross-language test vectors
    ├── constants/
    └── test_vectors/
```

## Documentation

- [C++ library](cpp/README.md) — build, test, and layout
- [Python package](python/README.md) — install, test, and optional Ray parallel sweeps

## Quick start

From the repository root:

```bash
make test          # C++ tests
make python-test   # Python tests
```

## Keeping C++ and Python in sync

- Physical constants and lookup tables live in `shared/constants/`.
- Cross-language test cases live in `shared/test_vectors/`.
- When changing numerical behavior, update both implementations and the shared fixtures in the same change.