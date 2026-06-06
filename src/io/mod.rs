//! VTK export (vti velocity fields, vtp STL surfaces).

pub mod vtk;
pub mod vtp;

pub use vtk::write_vti_velocity;
pub use vtp::write_vtp_stl_surface;
