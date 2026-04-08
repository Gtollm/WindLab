//! Bhatnagar–Gross–Krook (BGK) collision operator

use rayon::prelude::*;

use crate::grid::SoaDomain;
use crate::grid::cell::NodeType;
use crate::lattice::Q;
use crate::lattice::equilibrium::feq_populations;

#[inline]
pub fn omega_from_tau(tau: f64) -> f64 {
    1.0 / tau
}

#[inline]
pub fn nu_from_tau(tau: f64) -> f64 {
    crate::lattice::CS2 * (tau - 0.5)
}

pub fn collide_soa(domain: &mut SoaDomain, omega: f64) {
    let n = domain.ncells();
    let addr = domain as *mut SoaDomain as usize;

    (0..n).into_par_iter().for_each(move |j| {
        let d = unsafe { &mut *(addr as *mut SoaDomain) };
        match d.node_type[j] {
            NodeType::Solid => {}
            NodeType::Fluid | NodeType::Outlet(_) => {
                let rho = d.rho[j];
                let ux = d.ux[j];
                let uy = d.uy[j];
                let uz = d.uz[j];
                let mut feq = [0.0; Q];
                feq_populations(rho, ux, uy, uz, &mut feq);
                for (i, feq_val) in feq.iter().enumerate().take(Q) {
                    d.f[i][j] -= omega * (d.f[i][j] - *feq_val);
                }
            }
            NodeType::Inlet(u) => {
                let rho = d.rho[j];
                let mut feq = [0.0; Q];
                feq_populations(rho, u.x, u.y, u.z, &mut feq);
                for (i, feq_val) in feq.iter().enumerate().take(Q) {
                    d.f[i][j] -= omega * (d.f[i][j] - *feq_val);
                }
            }
        }
    });
}
