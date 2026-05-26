//! Solid-fill voxelization of triangle meshes using x-ray parity.
//!
//! For each (y, z) grid column, casts a ray parallel to x, finds all
//! intersections with the surface triangles, and marks cells between an
//! odd number of crossings as solid (standard even-odd winding rule).
//! This correctly fills the interior of any closed, manifold STL mesh
//! regardless of triangle density or body orientation.

use nalgebra::Vector3;

use super::stl::{Bounds, Tri};
use crate::grid::cell::NodeType;

pub fn voxelize_triangles(
    nx: usize,
    ny: usize,
    nz: usize,
    tris: &[Tri],
    world_bounds: &Bounds,
) -> Vec<NodeType> {
    let n = nx * ny * nz;
    let mut out = vec![NodeType::Fluid; n];

    if tris.is_empty() || nx < 2 || ny < 2 || nz < 2 {
        return out;
    }

    let wmin = world_bounds.min;
    let wmax = world_bounds.max;
    let extent = wmax - wmin;
    let dx = extent.x / (nx - 1) as f64;
    let dy = extent.y / (ny - 1) as f64;
    let dz = extent.z / (nz - 1) as f64;

    // Precompute per-triangle YZ bounding boxes for fast ray rejection.
    let yz_boxes: Vec<(f64, f64, f64, f64)> = tris
        .iter()
        .map(|t| {
            let ylo = t.a.y.min(t.b.y).min(t.c.y);
            let yhi = t.a.y.max(t.b.y).max(t.c.y);
            let zlo = t.a.z.min(t.b.z).min(t.c.z);
            let zhi = t.a.z.max(t.b.z).max(t.c.z);
            (ylo, yhi, zlo, zhi)
        })
        .collect();

    for iz in 0..nz {
        let z = wmin.z + iz as f64 * dz;
        for iy in 0..ny {
            let y = wmin.y + iy as f64 * dy;

            // Collect all x-coordinates where the (y,z) ray crosses a triangle.
            let mut xs: Vec<f64> = Vec::new();
            for (tri, &(ylo, yhi, zlo, zhi)) in tris.iter().zip(yz_boxes.iter()) {
                if y < ylo || y > yhi || z < zlo || z > zhi {
                    continue;
                }
                if let Some(x) = ray_x_tri(y, z, tri) {
                    xs.push(x);
                }
            }

            if xs.is_empty() {
                continue;
            }

            xs.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

            // Parity fill: a cell is solid if an odd number of crossings
            // lie strictly before it (even-odd winding rule).
            let mut ci = 0usize;
            let mut parity: usize = 0;
            for ix in 0..nx {
                let x = wmin.x + ix as f64 * dx;
                while ci < xs.len() && xs[ci] < x {
                    parity += 1;
                    ci += 1;
                }
                if parity % 2 == 1 {
                    out[crate::lattice::index(nx, ny, ix, iy, iz)] = NodeType::Solid;
                }
            }
        }
    }

    out
}

/// Returns the x-coordinate where the x-parallel ray at (y, z) pierces
/// `tri`, or `None` if the ray misses the triangle.
///
/// Uses barycentric coordinates in the YZ projection:
///   P = A + u·(C−A) + v·(B−A)  with  u,v ≥ 0  and  u+v ≤ 1
/// then interpolates x = A.x + u·(C.x−A.x) + v·(B.x−A.x).
fn ray_x_tri(y: f64, z: f64, tri: &Tri) -> Option<f64> {
    let a = tri.a;
    let b = tri.b;
    let c = tri.c;

    let v0y = c.y - a.y;
    let v0z = c.z - a.z;
    let v1y = b.y - a.y;
    let v1z = b.z - a.z;
    let v2y = y - a.y;
    let v2z = z - a.z;

    let denom = v0y * v1z - v0z * v1y;
    if denom.abs() < 1e-30 {
        return None; // triangle is degenerate or parallel to x-axis
    }

    let u = (v2y * v1z - v2z * v1y) / denom;
    let v = (v0y * v2z - v0z * v2y) / denom;

    if u >= 0.0 && v >= 0.0 && u + v <= 1.0 {
        Some(a.x + u * (c.x - a.x) + v * (b.x - a.x))
    } else {
        None
    }
}

pub fn world_to_lattice(
    p: Vector3<f64>,
    world_min: Vector3<f64>,
    world_max: Vector3<f64>,
    nx: usize,
    ny: usize,
    nz: usize,
) -> (usize, usize, usize) {
    use super::coords::{world_to_node_index_component, world_to_normalized};
    let n = world_to_normalized(p, world_min, world_max);
    (
        world_to_node_index_component(n.x, nx),
        world_to_node_index_component(n.y, ny),
        world_to_node_index_component(n.z, nz),
    )
}
