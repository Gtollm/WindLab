//! Macroscopic moments and body forcing

pub mod forcing;
pub mod macroscopic;

pub use forcing::apply_force_soa;
pub use macroscopic::update_macroscopic_soa;
