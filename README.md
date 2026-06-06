# WindLab

3D lattice Boltzmann CFD solver in Rust: **D3Q19** BGK collide-stream on a structured grid, **Rayon**-parallel macro steps, **Tokio** for async VTK I/O, **STL** triangle voxelization for solid geometry, **VTK** (`.vti` / `.vtp`) export, and optional live **Rerun** visualization.

## Build

```bash
cargo build --release
```

The release binary is `target/release/windlab`.

## Configuration

Copy and edit the example TOML:

```bash
cp config.example.toml my_case.toml
```

See [`config.example.toml`](config.example.toml) for grid size, relaxation time `tau`, inlet velocity, boundary tags, STL path, VTK output cadence, and run length.

## CLI

All commands are subcommands of `windlab`:

```bash
cargo run --release -- <subcommand> [flags]
# or after building:
./target/release/windlab <subcommand> [flags]
```

### `run`

Config-driven simulation; writes velocity fields under `io.output_dir` (default `output/`).

```bash
cargo run --release -- run --config config.example.toml
```

| Flag | Default | Description |
|------|---------|-------------|
| `--config <path>` | *(required)* | TOML config file |
| `--async-io` | enabled | Non-blocking VTK writes |
| `--no-progress` | `false` | Hide progress bar |
| `--rerun` | `false` | Open Rerun viewer and stream fields |
| `--slice-z <spec>` | `nz/2` | Z-planes for slice arrows (see below) |

Rerun logging during `run` follows `run.vtk_every` in the config (same steps as VTK output).

### `demo-channel`

Poiseuille-style channel on an `n x n x n` grid (no config file).

```bash
cargo run --release -- demo-channel --n 32 --steps 4000
```

| Flag | Default | Description |
|------|---------|-------------|
| `--n <size>` | `32` | Grid resolution (cube) |
| `--steps <N>` | `2000` | LBM steps |
| `--no-progress` | `false` | Hide progress bar |
| `--rerun` | `false` | Enable Rerun viewer |
| `--viz-every <N>` | `steps/20` | Rerun log interval (steps) |
| `--slice-z <spec>` | `nz/2` | Z-planes to visualize |

### `demo-sphere`

Flow past a sphere; compares simulated drag to empirical correlation.

```bash
cargo run --release -- demo-sphere --re 100 --diameter 20 --steps 3000
```

| Flag | Default | Description |
|------|---------|-------------|
| `--re <Re>` | `100` | Reynolds number |
| `--diameter <cells>` | `20` | Sphere diameter in cells |
| `--steps <N>` | `3000` | LBM steps |
| `--tau <tau>` | `0.6` | Relaxation time |
| `--no-progress` | `false` | Hide progress bar |
| `--rerun` | `false` | Enable Rerun viewer |

### `eval-stl`

Fast drag evaluation for an arbitrary STL (used by the shape optimizer). Builds a padded grid around the mesh, runs the solver, and prints a single line on **stdout**:

```
Cd=<value>
```

The optimizer parses lines that start with `Cd=`. Keep stdout clean aside from that contract.

```bash
cargo run --release -- eval-stl --stl path/to/body.stl --re 50 --steps 8000
```

| Flag | Default | Description |
|------|---------|-------------|
| `--stl <path>` | *(required)* | Input STL |
| `--re <Re>` | `50` | Reynolds number |
| `--tau <tau>` | `0.6` | Relaxation time |
| `--cpd <N>` | `30` | Cells per cross-stream diameter |
| `--steps <N>` | `8000` | LBM steps |
| `--no-progress` | `false` | Hide progress bar |
| `--quiet` | `false` | Suppress stderr diagnostics (optimizer-friendly) |

## Live visualization (Rerun)

Rerun is built in; pass `--rerun` on `run`, `demo-channel`, or `demo-sphere`.

Example channel demo with explicit planes and log interval:

```bash
cargo run --release -- demo-channel --n 64 --steps 800 --rerun --viz-every 10 --slice-z "16,24,28-31"
```

The viewer shows velocity arrows on selected Z-planes, a 3D velocity point cloud, solid geometry, and drag scalars over time.

### `--slice-z` syntax

| Input | Meaning |
|-------|---------|
| `16` | Single plane at z = 16 |
| `1,5,10` | Three planes |
| `3-8` | Planes z = 3 .. 8 (inclusive) |
| `1,2,3-8,10` | Mixed list and ranges |

If `--slice-z` is omitted, the middle plane `nz/2` is used.

### `--viz-every`

Only on `demo-channel`. Log to Rerun every **N** steps. Default: `max(1, steps/20)`.

For `run --rerun`, set `vtk_every` in the TOML config instead.

## Shape optimizer

Python **CMA-ES** search over axisymmetric bodies to minimize or maximize **Cd** at a given Reynolds number.

```bash
pip install -r optimizer/requirements.txt
cargo build --release

python optimizer/optimize.py --re 50 --mode min   # minimize drag
python optimizer/optimize.py --re 50 --mode max   # maximize drag
```

Each evaluation calls `windlab eval-stl` and reads `Cd=` from stdout. Iteration STLs are written under `runs/` (gitignored).

Optional flags: `--volume`, `--length`, `--max-iter`, etc. See `optimizer/optimize.py --help`.

For config-driven Rerun runs with geometry from TOML, use solid wall boundaries and a non-zero `inlet_u` where appropriate (see `optimizer/viz.toml` or `car.toml` in the repo).

## Physics notes

- Kinematic viscosity: `nu = (tau - 0.5) / 3` (lattice units).
- Inlet speed from Reynolds number: `u_inlet = Re * nu / cpd` (used inside `eval-stl`).

## License

MIT
