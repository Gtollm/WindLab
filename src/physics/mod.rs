//! Macroscopic moments, body forcing, and aerodynamic coefficients

pub mod drag;
pub mod forcing;
pub mod macroscopic;

pub use drag::{drag_coefficient, drag_force, projected_area_yz, reference_velocity};
pub use forcing::apply_force_soa;
pub use macroscopic::update_macroscopic_soa;
