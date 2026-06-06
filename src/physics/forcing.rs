//! Body force added after collision.

use rayon::prelude::*;

use crate::grid::NodeType;
use crate::grid::SoaDomain;
use crate::lattice::{C, CS2, Q, W};

pub fn apply_force_soa(domain: &mut SoaDomain, gx: f64, gy: f64, gz: f64, omega: f64) {
    let pref = 1.0 - 0.5 * omega;
    let n = domain.ncells();
    let addr = domain as *mut SoaDomain as usize;

    (0..n).into_par_iter().for_each(move |j| {
        let d = unsafe { &mut *(addr as *mut SoaDomain) };
        if matches!(d.node_type[j], NodeType::Solid) {
            return;
        }

        let ux = d.ux[j];
        let uy = d.uy[j];
        let uz = d.uz[j];
        let u_dot_g = ux * gx + uy * gy + uz * gz;

        for i in 0..Q {
            let ex = C[i][0] as f64;
            let ey = C[i][1] as f64;
            let ez = C[i][2] as f64;
            let eu = ex * ux + ey * uy + ez * uz;
            let eg = ex * gx + ey * gy + ez * gz;
            let fi = W[i] * (3.0 * eg + 9.0 * eu * eg / CS2 - 3.0 * u_dot_g);
            d.f[i][j] += pref * fi;
        }
    });
}
