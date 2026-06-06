//! STL mesh loading and conservative voxelization.

pub mod coords;
pub mod stl;
pub mod voxelizer;

pub use coords::world_to_normalized;
pub use stl::{load_stl_triangles, rotate_tris, Bounds, Tri};
pub use voxelizer::voxelize_triangles;
