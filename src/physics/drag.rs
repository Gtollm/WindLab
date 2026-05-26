use crate::grid::cell::NodeType;
use crate::grid::SoaDomain;
use crate::lattice::{C, OPPOSITE, Q};

pub fn drag_force(domain: &SoaDomain) -> [f64; 3] {
    let nx = domain.nx;
    let ny = domain.ny;
    let nz = domain.nz;

    let mut fx = 0.0_f64;
    let mut fy = 0.0_f64;
    let mut fz = 0.0_f64;

    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let id = domain.idx(x, y, z);
                if matches!(domain.node_type[id], NodeType::Solid) {
                    continue;
                }
                for i in 1..Q {
                    let nx_ = x as i32 + C[i][0];
                    let ny_ = y as i32 + C[i][1];
                    let nz_ = z as i32 + C[i][2];

                    if nx_ < 0 || ny_ < 0 || nz_ < 0
                        || nx_ >= nx as i32
                        || ny_ >= ny as i32
                        || nz_ >= nz as i32
                    {
                        continue;
                    }

                    let nb = domain.idx(nx_ as usize, ny_ as usize, nz_ as usize);
                    if matches!(domain.node_type[nb], NodeType::Solid) {
                        // Use the bounce-back population (OPPOSITE[i]) which equals
                        // the pre-streaming f[i] from the previous step.
                        // After pull streaming: f[OPPOSITE[i]][xf] = old f[i][xf].
                        let fi = domain.f[OPPOSITE[i]][id];
                        fx += 2.0 * fi * C[i][0] as f64;
                        fy += 2.0 * fi * C[i][1] as f64;
                        fz += 2.0 * fi * C[i][2] as f64;
                    }
                }
            }
        }
    }

    [fx, fy, fz]
}

#[inline]
pub fn drag_coefficient(f_drag: f64, rho_ref: f64, u_ref: f64, a_ref: f64) -> f64 {
    2.0 * f_drag / (rho_ref * u_ref * u_ref * a_ref)
}

pub fn projected_area_yz(domain: &SoaDomain) -> usize {
    let nx = domain.nx;
    let ny = domain.ny;
    let nz = domain.nz;
    let mut count = 0usize;
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if matches!(domain.node_type[domain.idx(x, y, z)], NodeType::Solid) {
                    count += 1;
                    break; // one hit per (y,z) column is enough
                }
            }
        }
    }
    count
}

pub fn reference_velocity(domain: &SoaDomain) -> f64 {
    let mut sum = 0.0_f64;
    let mut count = 0usize;
    for j in 0..domain.ncells() {
        match &domain.node_type[j] {
            NodeType::Inlet(u) => {
                sum += (u.x * u.x + u.y * u.y + u.z * u.z).sqrt();
                count += 1;
            }
            _ => {}
        }
    }
    if count > 0 {
        return sum / count as f64;
    }
    
    let mut sum = 0.0_f64;
    let mut n = 0usize;
    for j in 0..domain.ncells() {
        if !matches!(domain.node_type[j], NodeType::Solid) {
            sum += (domain.ux[j] * domain.ux[j]
                + domain.uy[j] * domain.uy[j]
                + domain.uz[j] * domain.uz[j])
                .sqrt();
            n += 1;
        }
    }
    if n > 0 { sum / n as f64 } else { 1.0 }
}
