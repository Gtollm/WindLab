WindLab

3D lattice Boltzmann CFD (**D3Q19**, BGK) in Rust: **Rayon**-parallel collide/stream/macro steps, **Tokio** for async
VTK export, optional **STL** voxelization, optional **Bevy** preview (`--features visualize`).

## Build

```bash
cargo build --release
```

Optional visualization (pulls Bevy; requires a compatible Rust toolchain):

```bash
cargo build --features visualize --release
```

## CLI

```bash
# Example config-driven run (writes `output/*.vti`)
cargo run --release -- run --config config.example.toml

# Quick channel demo
cargo run --release -- demo-channel --n 32 --steps 4000
```

With visualization enabled, set `WINDLAB_PREVIEW=1` to open a Bevy window while running from config (simulation advances
on a background thread).
