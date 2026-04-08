//! Structure-of-Arrays domain

use nalgebra::Vector3;

use super::cell::NodeType;
use super::in_bounds;
use crate::lattice::{equilibrium::feq_populations, index, Q};

#[derive(Clone, Debug)]
pub struct SoaDomain {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub dx: f64,
    pub dt: f64,
    pub f: [Vec<f64>; Q],
    pub f_tmp: [Vec<f64>; Q],
    pub node_type: Vec<NodeType>,
    pub rho: Vec<f64>,
    pub ux: Vec<f64>,
    pub uy: Vec<f64>,
    pub uz: Vec<f64>,
}

impl SoaDomain {
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64, dt: f64) -> Self {
        let n = nx * ny * nz;
        let mut eq = [0.0; Q];
        feq_populations(1.0, 0.0, 0.0, 0.0, &mut eq);

        Self {
            nx,
            ny,
            nz,
            dx,
            dt,
            f: std::array::from_fn(|_| vec![eq[0]; n]),
            f_tmp: std::array::from_fn(|_| vec![eq[0]; n]),
            node_type: vec![NodeType::Fluid; n],
            rho: vec![1.0; n],
            ux: vec![0.0; n],
            uy: vec![0.0; n],
            uz: vec![0.0; n],
        }
    }

    #[inline]
    pub fn ncells(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    #[inline]
    pub fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        index(self.nx, self.ny, x, y, z)
    }

    #[inline]
    pub fn in_bounds(&self, x: i32, y: i32, z: i32) -> bool {
        in_bounds(self.nx, self.ny, self.nz, x, y, z)
    }

    pub fn init_uniform(&mut self, rho0: f64, ux0: f64, uy0: f64, uz0: f64) {
        let n = self.ncells();
        let mut eq = [0.0; Q];
        feq_populations(rho0, ux0, uy0, uz0, &mut eq);
        for (i, eq_val) in eq.iter().enumerate().take(Q) {
            self.f[i].fill(*eq_val);
            self.f_tmp[i].fill(*eq_val);
        }
        for j in 0..n {
            if matches!(
                self.node_type[j],
                NodeType::Fluid | NodeType::Inlet(_) | NodeType::Outlet(_)
            ) {
                self.rho[j] = rho0;
                self.ux[j] = ux0;
                self.uy[j] = uy0;
                self.uz[j] = uz0;
            }
        }
    }

    pub fn copy_node_types_from(&mut self, types: &[NodeType]) {
        debug_assert_eq!(types.len(), self.node_type.len());
        self.node_type.copy_from_slice(types);
    }

    pub fn physical_velocity(&self, lattice_u: Vector3<f64>) -> Vector3<f64> {
        lattice_u * (self.dx / self.dt)
    }
}
