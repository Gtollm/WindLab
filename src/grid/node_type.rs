//! Per-node storage types.

use nalgebra::Vector3;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum NodeType {
    #[default]
    Fluid,
    Solid,
    Inlet(Vector3<f64>),
    Outlet(f64),
}
