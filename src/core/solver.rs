//! LBM time integration orchestrator

use crate::boundary::zou_he::apply_zou_he_soa;
use crate::collision::{collide_soa, omega_from_tau};
use crate::grid::SoaDomain;
use crate::physics::{apply_force_soa, update_macroscopic_soa};
use crate::streaming::stream_soa;

#[derive(Clone, Debug, Copy)]
pub struct LbmParams {
    pub tau: f64,
    pub body_force: [f64; 3],
}

impl Default for LbmParams {
    fn default() -> Self {
        Self {
            tau: 0.8,
            body_force: [0.0; 3],
        }
    }
}

pub fn step_soa(domain: &mut SoaDomain, p: &LbmParams) {
    let omega = omega_from_tau(p.tau);

    collide_soa(domain, omega);

    let g = p.body_force;
    if is_force_nonzero(g) {
        apply_force_soa(domain, g[0], g[1], g[2], omega);
    }

    stream_soa(domain);
    apply_zou_he_soa(domain);
    update_macroscopic_soa(domain);
}

pub fn run_soa(domain: &mut SoaDomain, p: &LbmParams, n_steps: usize) {
    for _ in 0..n_steps {
        step_soa(domain, p);
    }
}

#[inline]
pub fn nu_lattice(tau: f64) -> f64 {
    crate::collision::nu_from_tau(tau)
}

#[inline]
fn is_force_nonzero(g: [f64; 3]) -> bool {
    g[0].abs() > f64::EPSILON || g[1].abs() > f64::EPSILON || g[2].abs() > f64::EPSILON
}
