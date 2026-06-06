//! Macroscopic moment recovery: density and velocity from populations.

use rayon::prelude::*;

use crate::grid::NodeType;
use crate::grid::SoaDomain;
use crate::lattice::{C, Q};

pub fn update_macroscopic_soa(domain: &mut SoaDomain) {
    let n = domain.ncells();
    let addr = domain as *mut SoaDomain as usize;

    (0..n).into_par_iter().for_each(move |j| {
        let d = unsafe { &mut *(addr as *mut SoaDomain) };
        if matches!(d.node_type[j], NodeType::Solid) {
            d.rho[j] = 1.0;
            d.ux[j] = 0.0;
            d.uy[j] = 0.0;
            d.uz[j] = 0.0;
            return;
        }

        let mut rho = 0.0_f64;
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        let mut mz = 0.0_f64;

        for (i, c) in C.iter().enumerate().take(Q) {
            let fi = d.f[i][j];
            rho += fi;
            mx += c[0] as f64 * fi;
            my += c[1] as f64 * fi;
            mz += c[2] as f64 * fi;
        }

        let inv = 1.0 / rho.max(1e-12);
        d.rho[j] = rho;
        d.ux[j] = mx * inv;
        d.uy[j] = my * inv;
        d.uz[j] = mz * inv;
    });
}
