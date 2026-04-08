//! Link-wise bounce-back for no-slip walls

use crate::grid::cell::NodeType;
use crate::grid::Domain;
use crate::lattice::{index, unravel_index, C, OPPOSITE, Q};

pub fn pull_populations(domain: &Domain, id: usize, out: &mut [f64; Q]) {
    let nx = domain.nx;
    let ny = domain.ny;
    let nz = domain.nz;
    let (x, y, z) = unravel_index(nx, ny, id);
    let self_f = domain.cells[id].f;

    for i in 0..Q {
        let sx = x as i32 - C[i][0];
        let sy = y as i32 - C[i][1];
        let sz = z as i32 - C[i][2];

        out[i] = if sx >= 0
            && sy >= 0
            && sz >= 0
            && (sx as usize) < nx
            && (sy as usize) < ny
            && (sz as usize) < nz
        {
            let sid = index(nx, ny, sx as usize, sy as usize, sz as usize);
            if matches!(domain.cells[sid].node_type, NodeType::Solid) {
                self_f[OPPOSITE[i]]
            } else {
                domain.cells[sid].f[i]
            }
        } else {
            self_f[OPPOSITE[i]]
        };
    }
}

pub fn clear_solids(domain: &mut Domain) {
    for c in domain.cells.iter_mut() {
        if matches!(c.node_type, NodeType::Solid) {
            c.f.fill(0.0);
            c.f_tmp.fill(0.0);
        }
    }
}
