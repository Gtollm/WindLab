//! Structured 3D lattice domain with Structure-of-Arrays layout.

mod soa_domain;

pub mod block;
pub mod node_type;

pub use node_type::NodeType;
pub use soa_domain::SoaDomain;

#[inline]
pub fn in_bounds(nx: usize, ny: usize, nz: usize, x: i32, y: i32, z: i32) -> bool {
    x >= 0 && y >= 0 && z >= 0 && (x as usize) < nx && (y as usize) < ny && (z as usize) < nz
}
