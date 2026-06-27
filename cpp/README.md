# C++ library

Static library and tests for interpolation (`alpha`, `beta`, `gamma`) and banded Gaussian elimination.

## Layout

```
cpp/
├── include/    # headers
├── src/        # implementation
└── tests/      # test binaries and helpers
```

## Build and test

From this directory:

```bash
make          # builds libspadsolver.a
make test     # runs interpolation and gaussian elimination tests
make clean    # removes build artifacts
```

Or from the repository root:

```bash
make test
```

## Targets

| Command | Result |
|---------|--------|
| `make` | Build `libspadsolver.a` from `src/*.cpp` |
| `make test` | Build and run `test_interpolation` and `test_gaussian_elimination` |
| `make test-interpolation` | Interpolation tests only |
| `make test-gaussian-elimination` | Gaussian elimination tests only |
| `make clean` | Remove `.o`, `.a`, and test binaries |
