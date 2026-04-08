//! Structured 3D lattice domain with Structure-of-Arrays layout

mod soa_domain;

pub mod block;
pub mod cell;

pub use cell::{Cell, NodeType};
pub use soa_domain::SoaDomain;

#[inline]
pub fn in_bounds(nx: usize, ny: usize, nz: usize, x: i32, y: i32, z: i32) -> bool {
    x >= 0 && y >= 0 && z >= 0 && (x as usize) < nx && (y as usize) < ny && (z as usize) < nz
}

use crate::lattice::index;

#[derive(Clone, Debug)]
pub struct Domain {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub dx: f64,
    pub dt: f64,
    pub cells: Vec<Cell>,
}

impl Domain {
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64, dt: f64) -> Self {
        let n = nx * ny * nz;
        let cells = (0..n).map(|_| Cell::new_fluid_uniform(1.0)).collect();
        Self {
            nx,
            ny,
            nz,
            dx,
            dt,
            cells,
        }
    }

    #[inline]
    pub fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        index(self.nx, self.ny, x, y, z)
    }

    #[inline]
    pub fn in_bounds(&self, x: i32, y: i32, z: i32) -> bool {
        in_bounds(self.nx, self.ny, self.nz, x, y, z)
    }
}
