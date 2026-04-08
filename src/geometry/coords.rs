//! Coordinate transforms between world, normalised, and lattice index spaces

use nalgebra::Vector3;

#[inline]
pub fn nm1(n: usize) -> f64 {
    n.saturating_sub(1).max(1) as f64
}

#[inline]
pub fn norm_span_to_index_range(t_lo: f64, t_hi: f64, n: usize) -> (usize, usize) {
    let span = nm1(n);
    let (lo, hi) = if t_lo <= t_hi {
        (t_lo.clamp(0.0, 1.0), t_hi.clamp(0.0, 1.0))
    } else {
        (t_hi.clamp(0.0, 1.0), t_lo.clamp(0.0, 1.0))
    };

    let mut i0 = (lo * span).floor() as isize;
    let mut i1 = (hi * span).ceil() as isize;
    let imax = n.saturating_sub(1) as isize;

    i0 = i0.clamp(0, imax);
    i1 = i1.clamp(0, imax);

    if i1 < i0 {
        (i0 as usize, i0 as usize)
    } else {
        (i0 as usize, i1 as usize)
    }
}

#[inline]
pub fn world_to_normalized(
    p: Vector3<f64>,
    world_min: Vector3<f64>,
    world_max: Vector3<f64>,
) -> Vector3<f64> {
    let extent = world_max - world_min;
    let inv = Vector3::new(
        1.0 / extent.x.max(1e-30),
        1.0 / extent.y.max(1e-30),
        1.0 / extent.z.max(1e-30),
    );
    let rel = p - world_min;
    Vector3::new(
        (rel.x * inv.x).clamp(0.0, 1.0),
        (rel.y * inv.y).clamp(0.0, 1.0),
        (rel.z * inv.z).clamp(0.0, 1.0),
    )
}

#[inline]
pub fn world_to_node_index_component(t: f64, n: usize) -> usize {
    let i = (t.clamp(0.0, 1.0) * nm1(n)).round() as isize;
    i.clamp(0, n.saturating_sub(1) as isize) as usize
}
