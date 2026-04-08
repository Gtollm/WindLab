//! Per-node storage types

use nalgebra::Vector3;

use crate::lattice::{equilibrium::feq_populations, Q};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum NodeType {
    #[default]
    Fluid,
    Solid,
    Inlet(Vector3<f64>),
    Outlet(f64),
}

#[derive(Clone, Debug)]
pub struct Cell {
    pub f: [f64; Q],
    pub f_tmp: [f64; Q],
    pub node_type: NodeType,
    pub rho: f64,
    pub ux: f64,
    pub uy: f64,
    pub uz: f64,
}

impl Cell {
    pub fn new_fluid_uniform(rho: f64) -> Self {
        let mut f = [0.0; Q];
        feq_populations(rho, 0.0, 0.0, 0.0, &mut f);
        Self {
            f,
            f_tmp: f,
            node_type: NodeType::Fluid,
            rho,
            ux: 0.0,
            uy: 0.0,
            uz: 0.0,
        }
    }

    pub fn set_equilibrium(&mut self, rho: f64, ux: f64, uy: f64, uz: f64) {
        feq_populations(rho, ux, uy, uz, &mut self.f);
        self.f_tmp = self.f;
        self.rho = rho;
        self.ux = ux;
        self.uy = uy;
        self.uz = uz;
    }
}
