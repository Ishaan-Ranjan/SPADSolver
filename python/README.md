# Python package

Python port of the SPADSolver numerical core: `alpha`, `beta`, `gamma` interpolation and Gaussian elimination helpers.

## Virtual environment

Use a venv to keep dependencies isolated from system Python.

### Create and activate

```bash
cd python
python3 -m venv .venv
```

**macOS / Linux:**

```bash
source .venv/bin/activate
```

**Windows (Command Prompt):**

```cmd
.venv\Scripts\activate.bat
```

**Windows (PowerShell):**

```powershell
.venv\Scripts\Activate.ps1
```

Your shell prompt should show `(.venv)` when the environment is active.

### Install dependencies

```bash
python3 -m pip install --upgrade pip
pip install -r requirements.txt
```

This performs an editable install of `spadsolver` and installs pytest. Dependencies are also declared in `pyproject.toml`.

For optional Ray parallel support:

```bash
pip install -r requirements-ray.txt
```

### Deactivate

When you are done working in the project:

```bash
deactivate
```

### Notes

- `.venv/` is gitignored — each clone creates its own environment.
- After pulling code changes, reinstall if dependencies change: `pip install -r requirements.txt`
- Alternative without requirements files: `pip install -e ".[dev]"` or `pip install -e ".[dev,ray]"`

## Run tests

With the venv activated:

```bash
pytest
```

Or from the repository root:

```bash
make python-test
```

## Optional: parallel sweeps with Ray

Requires `pip install -r requirements-ray.txt` (or `pip install -e ".[dev,ray]"`).

```bash
python3 examples/parameter_sweep.py
```

Use the Ray helpers in application code:

```python
from spadsolver.ray_tasks import parallel_parameter_sweep

points = [(0.1, 3e6, 298.0, 3.0), (0.2, 3e6, 298.0, 3.5)]
results = parallel_parameter_sweep(points, chunk_size=64)
```

Ray is optional — the core `spadsolver` package runs without it. See `spadsolver/ray_tasks.py` for `parallel_parameter_sweep`, `parallel_solve_systems`, and related helpers.

## Layout

```
python/
├── spadsolver/          # package source
├── tests/               # pytest suite
├── examples/            # Ray parameter sweep example
├── requirements.txt     # dev install (editable + pytest)
└── requirements-ray.txt # optional Ray dependency
```
