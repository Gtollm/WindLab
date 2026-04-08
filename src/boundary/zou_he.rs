//! Simplified Zou–He style velocity / pressure boundaries
use nalgebra::Vector3;

use crate::grid::cell::NodeType;
use crate::grid::SoaDomain;
use crate::lattice::equilibrium::feq_populations;
use crate::lattice::{index, unravel_index, Q};

pub fn apply_zou_he_soa(domain: &mut SoaDomain) {
    let nx = domain.nx;
    let ny = domain.ny;
    let _nz = domain.nz;
    let n = domain.ncells();
    let mut buf = [0.0; Q];

    for j in 0..n {
        match domain.node_type[j] {
            NodeType::Inlet(u) => {
                let rho = domain.rho[j].max(1e-6);
                feq_populations(rho, u.x, u.y, u.z, &mut buf);
                for (i, val) in buf.iter().enumerate().take(Q) {
                    domain.f[i][j] = *val;
                }
                domain.ux[j] = u.x;
                domain.uy[j] = u.y;
                domain.uz[j] = u.z;
            }
            NodeType::Outlet(rho_t) => {
                let (x, y, z) = unravel_index(nx, ny, j);
                let mut ux = domain.ux[j];
                let mut uy = domain.uy[j];
                let mut uz = domain.uz[j];

                if x > 0 {
                    let nj = index(nx, ny, x - 1, y, z);
                    if matches!(
                        domain.node_type[nj],
                        NodeType::Fluid | NodeType::Inlet(_) | NodeType::Outlet(_)
                    ) {
                        ux = domain.ux[nj];
                        uy = domain.uy[nj];
                        uz = domain.uz[nj];
                    }
                }

                feq_populations(rho_t, ux, uy, uz, &mut buf);
                for (i, val) in buf.iter().enumerate().take(Q) {
                    domain.f[i][j] = *val;
                }
                domain.rho[j] = rho_t;
                domain.ux[j] = ux;
                domain.uy[j] = uy;
                domain.uz[j] = uz;
            }
            _ => {}
        }
    }
}

pub fn tag_x0_inlet(types: &mut [NodeType], nx: usize, ny: usize, nz: usize, u: Vector3<f64>) {
    for z in 0..nz {
        for y in 0..ny {
            let id = index(nx, ny, 0, y, z);
            if matches!(types[id], NodeType::Fluid) {
                types[id] = NodeType::Inlet(u);
            }
        }
    }
}

pub fn tag_xmax_outlet(types: &mut [NodeType], nx: usize, ny: usize, nz: usize, rho: f64) {
    let xmax = nx - 1;
    for z in 0..nz {
        for y in 0..ny {
            let id = index(nx, ny, xmax, y, z);
            if matches!(types[id], NodeType::Fluid) {
                types[id] = NodeType::Outlet(rho);
            }
        }
    }
}
