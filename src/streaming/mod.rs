//! Pull-based streaming with link-wise bounce-back at solid boundaries

use rayon::prelude::*;

use crate::grid::cell::NodeType;
use crate::grid::SoaDomain;
use crate::lattice::{index, unravel_index, C, OPPOSITE, Q};

pub fn stream_soa(domain: &mut SoaDomain) {
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

        let (x, y, z) = unravel_index(nx, ny, id);

        for i in 0..Q {
            let sx = x as i32 - C[i][0];
            let sy = y as i32 - C[i][1];
            let sz = z as i32 - C[i][2];

            d.f_tmp[i][id] = if sx >= 0
                && sy >= 0
                && sz >= 0
                && (sx as usize) < nx
                && (sy as usize) < ny
                && (sz as usize) < nz
            {
                let sid = index(nx, ny, sx as usize, sy as usize, sz as usize);
                if matches!(d.node_type[sid], NodeType::Solid) {
                    d.f[OPPOSITE[i]][id]
                } else {
                    d.f[i][sid]
                }
            } else {
                d.f[OPPOSITE[i]][id]
            };
        }
    });

    for i in 0..Q {
        std::mem::swap(&mut domain.f[i], &mut domain.f_tmp[i]);
    }
}
