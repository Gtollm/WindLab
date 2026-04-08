//! Conservative AABB-based voxelization of triangle meshes.

use nalgebra::Vector3;

use super::coords::{norm_span_to_index_range, world_to_node_index_component, world_to_normalized};
use super::stl::{Bounds, Tri};
use crate::grid::cell::NodeType;

fn tri_aabb(t: &Tri) -> (Vector3<f64>, Vector3<f64>) {
    let min = Vector3::new(
        t.a.x.min(t.b.x).min(t.c.x),
        t.a.y.min(t.b.y).min(t.c.y),
        t.a.z.min(t.b.z).min(t.c.z),
    );
    let max = Vector3::new(
        t.a.x.max(t.b.x).max(t.c.x),
        t.a.y.max(t.b.y).max(t.c.y),
        t.a.z.max(t.b.z).max(t.c.z),
    );
    (min, max)
}

pub fn voxelize_triangles(
    nx: usize,
    ny: usize,
    nz: usize,
    tris: &[Tri],
    world_bounds: &Bounds,
) -> Vec<NodeType> {
    let n = nx * ny * nz;
    let mut out = vec![NodeType::Fluid; n];
    let wmin = world_bounds.min;
    let wmax = world_bounds.max;

    for t in tris {
        let (tmin, tmax) = tri_aabb(t);
        let n0 = world_to_normalized(tmin, wmin, wmax);
        let n1 = world_to_normalized(tmax, wmin, wmax);

        let (xa, xb) = norm_span_to_index_range(n0.x, n1.x, nx);
        let (ya, yb) = norm_span_to_index_range(n0.y, n1.y, ny);
        let (za, zb) = norm_span_to_index_range(n0.z, n1.z, nz);

        for iz in za..=zb {
            for iy in ya..=yb {
                for ix in xa..=xb {
                    let id = crate::lattice::index(nx, ny, ix, iy, iz);
                    out[id] = NodeType::Solid;
                }
            }
        }
    }

    out
}

pub fn world_to_lattice(
    p: Vector3<f64>,
    world_min: Vector3<f64>,
    world_max: Vector3<f64>,
    nx: usize,
    ny: usize,
    nz: usize,
) -> (usize, usize, usize) {
    let n = world_to_normalized(p, world_min, world_max);
    (
        world_to_node_index_component(n.x, nx),
        world_to_node_index_component(n.y, ny),
        world_to_node_index_component(n.z, nz),
    )
}
