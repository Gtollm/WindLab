//! Body force - uniform body force density added after collision.

use crate::grid::cell::NodeType;
use crate::grid::SoaDomain;
use crate::lattice::{C, CS2, Q, W};

pub fn apply_force_soa(domain: &mut SoaDomain, gx: f64, gy: f64, gz: f64, omega: f64) {
    let pref = 1.0 - 0.5 * omega;
    let n = domain.ncells();

    for j in 0..n {
        if matches!(domain.node_type[j], NodeType::Solid) {
            continue;
        }

        let ux = domain.ux[j];
        let uy = domain.uy[j];
        let uz = domain.uz[j];
        let u_dot_g = ux * gx + uy * gy + uz * gz;

        for i in 0..Q {
            let ex = C[i][0] as f64;
            let ey = C[i][1] as f64;
            let ez = C[i][2] as f64;
            let eu = ex * ux + ey * uy + ez * uz;
            let eg = ex * gx + ey * gy + ez * gz;
            let fi = W[i] * (3.0 * eg + 9.0 * eu * eg / CS2 - 3.0 * u_dot_g);
            domain.f[i][j] += pref * fi;
        }
    }
}
