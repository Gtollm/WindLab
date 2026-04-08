//! Equilibrium distribution for D3Q19

use super::{C, CS2, Q, W};

#[inline(always)]
pub fn feq_populations(rho: f64, ux: f64, uy: f64, uz: f64, out: &mut [f64; Q]) {
    let u2 = ux * ux + uy * uy + uz * uz;
    let inv_cs2 = 1.0 / CS2;
    let u2_term = u2 * inv_cs2 * 0.5;

    for i in 0..Q {
        let eu = C[i][0] as f64 * ux + C[i][1] as f64 * uy + C[i][2] as f64 * uz;
        let eu_cs2 = eu * inv_cs2;
        out[i] = W[i] * rho * (1.0 + eu_cs2 + 0.5 * eu * eu_cs2 * inv_cs2 - u2_term);
    }
}

#[inline]
pub fn feq_i(i: usize, rho: f64, ux: f64, uy: f64, uz: f64) -> f64 {
    let eu = C[i][0] as f64 * ux + C[i][1] as f64 * uy + C[i][2] as f64 * uz;
    let u2 = ux * ux + uy * uy + uz * uz;
    let inv_cs2 = 1.0 / CS2;
    W[i] * rho * (1.0 + eu * inv_cs2 + 0.5 * eu * eu * inv_cs2 * inv_cs2 - u2 * inv_cs2 * 0.5)
}
