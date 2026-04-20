## Build

```bash
cargo build --release
```

Rerun visualization included by default — no extra feature flags needed.
## Live Visualization (Rerun)

```bash
cargo run --release -- demo-channel \
  --n 64 --steps 2000 --rerun --viz-every 10 --slice-z "16,24,28-31"
```

Opens the Rerun viewer automatically. Shows:
- 3D velocity arrows on selected Z-planes
- Full 3D velocity point cloud
- Solid geometry
- Convergence scalar (mean speed) over time

### `--slice-z` syntax

| Input | Meaning |
|-------|---------|
| `16` | single plane at z=16 |
| `1,5,10` | three planes |
| `3-8` | planes z=3,4,5,6,7,8 |
| `1,2,3-8,10` | mixed |

Omit `--slice-z` → defaults to `nz/2`.

### `--viz-every N`

Log to Rerun every N steps. Default: `steps/20`.

---

## CLI Reference

### `run`

| Flag | Default | Description |
|------|---------|-------------|
| `--config <path>` | required | TOML config file |
| `--async-io` | true | Non-blocking VTK writes |
| `--no-progress` | false | Hide progress bar |
| `--rerun` | false | Enable live Rerun viewer |
| `--slice-z <spec>` | nz/2 | Z-planes to visualize |

### `demo-channel`

| Flag | Default | Description |
|------|---------|-------------|
| `--n <size>` | 32 | Grid size (n×n×n) |
| `--steps <N>` | 2000 | Simulation steps |
| `--no-progress` | false | Hide progress bar |
| `--rerun` | false | Enable live Rerun viewer |
| `--viz-every <N>` | steps/20 | Rerun log interval |
| `--slice-z <spec>` | nz/2 | Z-planes to visualize |