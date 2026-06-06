//! Domain outer-face rules for pull streaming (solid / open / periodic).

use rayon::prelude::*;

use crate::config::{BoundaryConfig, BoundaryKind};
use crate::grid::NodeType;
use crate::grid::SoaDomain;
use crate::lattice::{index, C, OPPOSITE, Q};

use super::periodic::wrap;

#[inline]
fn in_domain(gx: i32, gy: i32, gz: i32, nx: usize, ny: usize, nz: usize) -> bool {
    gx >= 0 && gy >= 0 && gz >= 0 && (gx as usize) < nx && (gy as usize) < ny && (gz as usize) < nz
}

fn apply_periodic(
    mut gx: i32,
    mut gy: i32,
    mut gz: i32,
    nx: usize,
    ny: usize,
    nz: usize,
    b: &BoundaryConfig,
) -> (i32, i32, i32) {
    for _ in 0..8 {
        let mut changed = false;
        if gx < 0 && b.left == BoundaryKind::Periodic {
            gx = wrap(gx, nx) as i32;
            changed = true;
        }
        if gx >= nx as i32 && b.right == BoundaryKind::Periodic {
            gx = wrap(gx, nx) as i32;
            changed = true;
        }
        if gy < 0 && b.front == BoundaryKind::Periodic {
            gy = wrap(gy, ny) as i32;
            changed = true;
        }
        if gy >= ny as i32 && b.back == BoundaryKind::Periodic {
            gy = wrap(gy, ny) as i32;
            changed = true;
        }
        if gz < 0 && b.bottom == BoundaryKind::Periodic {
            gz = wrap(gz, nz) as i32;
            changed = true;
        }
        if gz >= nz as i32 && b.top == BoundaryKind::Periodic {
            gz = wrap(gz, nz) as i32;
            changed = true;
        }
        if !changed {
            break;
        }
    }
    (gx, gy, gz)
}

fn all_violated_faces_open(
    gx: i32,
    gy: i32,
    gz: i32,
    nx: usize,
    ny: usize,
    nz: usize,
    b: &BoundaryConfig,
) -> bool {
    let mut any = false;
    if gx < 0 {
        any = true;
        if b.left != BoundaryKind::Open {
            return false;
        }
    }
    if gx >= nx as i32 {
        any = true;
        if b.right != BoundaryKind::Open {
            return false;
        }
    }
    if gy < 0 {
        any = true;
        if b.front != BoundaryKind::Open {
            return false;
        }
    }
    if gy >= ny as i32 {
        any = true;
        if b.back != BoundaryKind::Open {
            return false;
        }
    }
    if gz < 0 {
        any = true;
        if b.bottom != BoundaryKind::Open {
            return false;
        }
    }
    if gz >= nz as i32 {
        any = true;
        if b.top != BoundaryKind::Open {
            return false;
        }
    }
    any
}

pub fn stream_soa_with_boundaries(domain: &mut SoaDomain, bounds: &BoundaryConfig) {
    let nx = domain.nx;
    let ny = domain.ny;
    let nz = domain.nz;
    let n = domain.ncells();
    let addr = domain as *mut SoaDomain as usize;

    (0..n).into_par_iter().for_each(move |id| {
        let d = unsafe { &mut *(addr as *mut SoaDomain) };
        if matches!(d.node_type[id], NodeType::Solid) {
            return;
        }

        let (x, y, z) = crate::lattice::unravel_index(nx, ny, id);
        let self_f: [f64; Q] = std::array::from_fn(|k| d.f[k][id]);

        for i in 0..Q {
            let sx = x as i32 - C[i][0];
            let sy = y as i32 - C[i][1];
            let sz = z as i32 - C[i][2];

            let (gx, gy, gz) = apply_periodic(sx, sy, sz, nx, ny, nz, bounds);

            d.f_tmp[i][id] = if in_domain(gx, gy, gz, nx, ny, nz) {
                let sid = index(nx, ny, gx as usize, gy as usize, gz as usize);
                if matches!(d.node_type[sid], NodeType::Solid) {
                    self_f[OPPOSITE[i]]
                } else {
                    d.f[i][sid]
                }
            } else if all_violated_faces_open(gx, gy, gz, nx, ny, nz, bounds) {
                self_f[i]
            } else {
                self_f[OPPOSITE[i]]
            };
        }
    });

    for i in 0..Q {
        std::mem::swap(&mut domain.f[i], &mut domain.f_tmp[i]);
    }
}

#[cfg(test)]
mod tests {
    use crate::config::{BoundaryConfig, BoundaryKind};

    use super::apply_periodic;

    #[test]
    fn periodic_wraps_negative_x() {
        let b = BoundaryConfig {
            left: BoundaryKind::Periodic,
            right: BoundaryKind::Solid,
            front: BoundaryKind::Solid,
            back: BoundaryKind::Solid,
            bottom: BoundaryKind::Solid,
            top: BoundaryKind::Solid,
        };
        let (gx, _, _) = apply_periodic(-1, 4, 4, 4, 4, 4, &b);
        assert_eq!(gx, 3);
    }
}
