//! Macroscopic moment recovery - density and velocity from populations

use crate::grid::cell::NodeType;
use crate::grid::SoaDomain;
use crate::lattice::{C, Q};

pub fn update_macroscopic_soa(domain: &mut SoaDomain) {
    let n = domain.ncells();

    for j in 0..n {
        if matches!(domain.node_type[j], NodeType::Solid) {
            domain.rho[j] = 1.0;
            domain.ux[j] = 0.0;
            domain.uy[j] = 0.0;
            domain.uz[j] = 0.0;
            continue;
        }

        let mut rho = 0.0_f64;
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        let mut mz = 0.0_f64;

        for (i, c) in C.iter().enumerate().take(Q) {
            let fi = domain.f[i][j];
            rho += fi;
            mx += c[0] as f64 * fi;
            my += c[1] as f64 * fi;
            mz += c[2] as f64 * fi;
        }

        let inv = 1.0 / rho.max(1e-12);
        domain.rho[j] = rho;
        domain.ux[j] = mx * inv;
        domain.uy[j] = my * inv;
        domain.uz[j] = mz * inv;
    }
}
