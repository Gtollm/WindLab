//! STL mesh loading and conservative voxelization

pub mod coords;
pub mod stl;
pub mod voxelizer;

pub use coords::{norm_span_to_index_range, world_to_normalized};
pub use stl::{load_stl_triangles, Bounds, Tri};
pub use voxelizer::{voxelize_triangles, world_to_lattice};
