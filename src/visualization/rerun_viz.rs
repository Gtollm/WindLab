//! Rerun visualization helpers for SoaDomain fields.
//!
//! Algorithms are untouched — this module only reads domain data and logs it.

use crate::geometry::stl::Tri;
use crate::grid::{cell::NodeType, SoaDomain};

pub use rerun::RecordingStream;

/// Spawn the Rerun viewer and return a recording stream.
/// Requires `rerun` CLI in PATH — install with: `cargo install rerun-cli --locked`
pub fn spawn_viewer(app_id: &str) -> Result<RecordingStream, rerun::RecordingStreamError> {
    rerun::RecordingStreamBuilder::new(app_id).spawn()
}

/// Log original STL triangles as a surface mesh — static, always visible.
///
/// Vertices are transformed from world space into grid index space using a
/// **uniform** scale (smallest axis wins, preserves aspect ratio) then centered
/// in the [0..nx-1] × [0..ny-1] × [0..nz-1] domain so they align with the
/// velocity field logged by `log_velocity_points`.
pub fn log_stl_mesh(
    rec: &RecordingStream,
    tris: &[Tri],
    bounds: &crate::geometry::stl::Bounds,
    nx: usize,
    ny: usize,
    nz: usize,
) -> Result<(), rerun::RecordingStreamError> {
    if tris.is_empty() {
        return Ok(());
    }

    let extent = bounds.max - bounds.min;
    // grid spans [0 .. n-1]; pick the tightest axis for uniform scale
    let scale = ((nx as f64 - 1.0) / extent.x.max(1e-30))
        .min((ny as f64 - 1.0) / extent.y.max(1e-30))
        .min((nz as f64 - 1.0) / extent.z.max(1e-30));

    // center the scaled model in the domain
    let ox = ((nx as f64 - 1.0) - extent.x * scale) * 0.5;
    let oy = ((ny as f64 - 1.0) - extent.y * scale) * 0.5;
    let oz = ((nz as f64 - 1.0) - extent.z * scale) * 0.5;

    let to_grid = |p: nalgebra::Vector3<f64>| -> [f32; 3] {
        [
            ((p.x - bounds.min.x) * scale + ox) as f32,
            ((p.y - bounds.min.y) * scale + oy) as f32,
            ((p.z - bounds.min.z) * scale + oz) as f32,
        ]
    };

    let verts: Vec<[f32; 3]> = tris
        .iter()
        .flat_map(|t| [to_grid(t.a), to_grid(t.b), to_grid(t.c)])
        .collect();

    let tri_indices: Vec<[u32; 3]> = (0..tris.len() as u32)
        .map(|i| [i * 3, i * 3 + 1, i * 3 + 2])
        .collect();

    let gray = rerun::Color::from_rgb(160, 160, 170);
    let colors = vec![gray; verts.len()];

    rec.log_static(
        "geometry/solid",
        &rerun::Mesh3D::new(verts)
            .with_triangle_indices(tri_indices)
            .with_vertex_colors(colors),
    )?;
    Ok(())
}

/// Log solid geometry as a surface mesh (only exposed faces) — static, always visible.
pub fn log_geometry(rec: &RecordingStream, domain: &SoaDomain) -> Result<(), rerun::RecordingStreamError> {
    let (verts, tris) = build_solid_surface(domain);
    if tris.is_empty() {
        return Ok(());
    }

    let gray = rerun::Color::from_rgb(130, 130, 140);
    let colors = vec![gray; verts.len()];

    rec.log_static(
        "geometry/solid",
        &rerun::Mesh3D::new(verts)
            .with_triangle_indices(tris)
            .with_vertex_colors(colors),
    )?;
    Ok(())
}

/// Extract only the faces of solid voxels that border a non-solid (or out-of-bounds) cell.
fn build_solid_surface(domain: &SoaDomain) -> (Vec<[f32; 3]>, Vec<[u32; 3]>) {
    let mut verts: Vec<[f32; 3]> = Vec::new();
    let mut tris: Vec<[u32; 3]> = Vec::new();

    let is_solid = |x: i32, y: i32, z: i32| -> bool {
        if x < 0
            || y < 0
            || z < 0
            || x >= domain.nx as i32
            || y >= domain.ny as i32
            || z >= domain.nz as i32
        {
            return false;
        }
        matches!(
            domain.node_type[domain.idx(x as usize, y as usize, z as usize)],
            NodeType::Solid
        )
    };

    for z in 0..domain.nz {
        for y in 0..domain.ny {
            for x in 0..domain.nx {
                let (ix, iy, iz) = (x as i32, y as i32, z as i32);
                if !is_solid(ix, iy, iz) {
                    continue;
                }
                let (fx, fy, fz) = (x as f32, y as f32, z as f32);

                // +X
                if !is_solid(ix + 1, iy, iz) {
                    add_quad(&mut verts, &mut tris,
                        [fx+1.0, fy,     fz    ],
                        [fx+1.0, fy,     fz+1.0],
                        [fx+1.0, fy+1.0, fz+1.0],
                        [fx+1.0, fy+1.0, fz    ]);
                }
                // -X
                if !is_solid(ix - 1, iy, iz) {
                    add_quad(&mut verts, &mut tris,
                        [fx, fy,     fz    ],
                        [fx, fy+1.0, fz    ],
                        [fx, fy+1.0, fz+1.0],
                        [fx, fy,     fz+1.0]);
                }
                // +Y
                if !is_solid(ix, iy + 1, iz) {
                    add_quad(&mut verts, &mut tris,
                        [fx,     fy+1.0, fz    ],
                        [fx+1.0, fy+1.0, fz    ],
                        [fx+1.0, fy+1.0, fz+1.0],
                        [fx,     fy+1.0, fz+1.0]);
                }
                // -Y
                if !is_solid(ix, iy - 1, iz) {
                    add_quad(&mut verts, &mut tris,
                        [fx,     fy, fz    ],
                        [fx,     fy, fz+1.0],
                        [fx+1.0, fy, fz+1.0],
                        [fx+1.0, fy, fz    ]);
                }
                // +Z
                if !is_solid(ix, iy, iz + 1) {
                    add_quad(&mut verts, &mut tris,
                        [fx,     fy,     fz+1.0],
                        [fx+1.0, fy,     fz+1.0],
                        [fx+1.0, fy+1.0, fz+1.0],
                        [fx,     fy+1.0, fz+1.0]);
                }
                // -Z
                if !is_solid(ix, iy, iz - 1) {
                    add_quad(&mut verts, &mut tris,
                        [fx,     fy,     fz    ],
                        [fx,     fy+1.0, fz    ],
                        [fx+1.0, fy+1.0, fz    ],
                        [fx+1.0, fy,     fz    ]);
                }
            }
        }
    }

    (verts, tris)
}

fn add_quad(
    verts: &mut Vec<[f32; 3]>,
    tris: &mut Vec<[u32; 3]>,
    v0: [f32; 3], v1: [f32; 3], v2: [f32; 3], v3: [f32; 3],
) {
    let base = verts.len() as u32;
    verts.extend_from_slice(&[v0, v1, v2, v3]);
    tris.push([base, base + 1, base + 2]);
    tris.push([base, base + 2, base + 3]);
}

/// Log fluid velocity as 3-D points colored by speed (turbo palette).
pub fn log_velocity_points(
    rec: &RecordingStream,
    domain: &SoaDomain,
    step: usize,
) -> Result<(), rerun::RecordingStreamError> {
    rec.set_time_sequence("step", step as i64);

    let v_max = max_speed(domain).max(1e-12);

    let mut positions: Vec<[f32; 3]> = Vec::with_capacity(domain.ncells());
    let mut colors: Vec<rerun::Color> = Vec::with_capacity(domain.ncells());

    for z in 0..domain.nz {
        for y in 0..domain.ny {
            for x in 0..domain.nx {
                let i = domain.idx(x, y, z);
                if is_fluid(&domain.node_type[i]) {
                    let speed = speed_at(domain, i);
                    positions.push([x as f32, y as f32, z as f32]);
                    colors.push(turbo((speed / v_max).clamp(0.0, 1.0) as f32));
                }
            }
        }
    }

    rec.log(
        "fluid/velocity",
        &rerun::Points3D::new(positions)
            .with_colors(colors)
            .with_radii([rerun::Radius::new_ui_points(2.5)]),
    )?;
    Ok(())
}

/// Log velocity arrows on one or more Z-slices (2-D cross sections, scaled for readability).
/// Log velocity arrows for the given Z-index planes.
/// `z_indices`: slice of z-coordinates to render (already clamped/resolved by caller).
pub fn log_velocity_slice(
    rec: &RecordingStream,
    domain: &SoaDomain,
    step: usize,
    z_indices: &[usize],
) -> Result<(), rerun::RecordingStreamError> {
    rec.set_time_sequence("step", step as i64);

    let v_max = max_speed(domain).max(1e-12);
    let arrow_scale = 0.4_f32 * (domain.nx.min(domain.ny) as f32) / v_max as f32;

    let mut origins: Vec<[f32; 3]> = Vec::new();
    let mut vectors: Vec<[f32; 3]> = Vec::new();
    let mut colors: Vec<rerun::Color> = Vec::new();

    for &z in z_indices {
        for y in 0..domain.ny {
            for x in 0..domain.nx {
                let i = domain.idx(x, y, z);
                if is_fluid(&domain.node_type[i]) {
                    let ux = domain.ux[i] as f32;
                    let uy = domain.uy[i] as f32;
                    let uz = domain.uz[i] as f32;
                    let speed = (domain.ux[i] * domain.ux[i]
                        + domain.uy[i] * domain.uy[i]
                        + domain.uz[i] * domain.uz[i])
                        .sqrt();
                    origins.push([x as f32, y as f32, z as f32]);
                    vectors.push([ux * arrow_scale, uy * arrow_scale, uz * arrow_scale]);
                    colors.push(turbo((speed / v_max).clamp(0.0, 1.0) as f32));
                }
            }
        }
    }

    if !origins.is_empty() {
        rec.log(
            "fluid/velocity_slice",
            &rerun::Arrows3D::from_vectors(vectors)
                .with_origins(origins)
                .with_colors(colors),
        )?;
    }
    Ok(())
}

// ── helpers ──────────────────────────────────────────────────────────────────

#[inline]
fn is_fluid(nt: &NodeType) -> bool {
    !matches!(nt, NodeType::Solid)
}

#[inline]
fn speed_at(domain: &SoaDomain, i: usize) -> f64 {
    let ux = domain.ux[i];
    let uy = domain.uy[i];
    let uz = domain.uz[i];
    (ux * ux + uy * uy + uz * uz).sqrt()
}

fn max_speed(domain: &SoaDomain) -> f64 {
    domain
        .ux
        .iter()
        .zip(domain.uy.iter())
        .zip(domain.uz.iter())
        .map(|((ux, uy), uz)| (ux * ux + uy * uy + uz * uz).sqrt())
        .fold(0.0_f64, f64::max)
}

/// Turbo color map: blue (0) → cyan → green → yellow → red (1).
fn turbo(t: f32) -> rerun::Color {
    let (r, g, b) = if t < 0.25 {
        let s = t * 4.0;
        (0u8, (s * 255.0) as u8, 255u8)
    } else if t < 0.5 {
        let s = (t - 0.25) * 4.0;
        (0, 255, ((1.0 - s) * 255.0) as u8)
    } else if t < 0.75 {
        let s = (t - 0.5) * 4.0;
        ((s * 255.0) as u8, 255, 0)
    } else {
        let s = (t - 0.75) * 4.0;
        (255, ((1.0 - s) * 255.0) as u8, 0)
    };
    rerun::Color::from_rgb(r, g, b)
}
