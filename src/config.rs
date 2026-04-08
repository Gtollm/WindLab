//! TOML configuration for simulations

use serde::Deserialize;
use std::path::Path;

use crate::core::solver::LbmParams;
use crate::grid::cell::NodeType;
use crate::grid::SoaDomain;

#[derive(Debug, Deserialize, Clone)]
pub struct SimConfig {
    pub grid: GridConfig,
    pub physics: PhysicsConfig,
    pub run: RunConfig,
    #[serde(default)]
    pub io: IoConfig,
    #[serde(default)]
    pub geometry: Option<GeometryConfig>,
}

#[derive(Debug, Deserialize, Clone)]
pub struct GridConfig {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub dx: f64,
    pub dt: f64,
}

#[derive(Debug, Deserialize, Clone)]
pub struct PhysicsConfig {
    pub tau: f64,
    #[serde(default)]
    pub body_force: [f64; 3],
}

impl PhysicsConfig {
    pub fn resolved_body_force(&self) -> [f64; 3] {
        self.body_force
    }
}

#[derive(Debug, Deserialize, Clone)]
pub struct RunConfig {
    pub steps: usize,
    #[serde(default = "default_vtk_every")]
    pub vtk_every: usize,
}

fn default_vtk_every() -> usize {
    100
}

#[derive(Debug, Deserialize, Clone, Default)]
pub struct IoConfig {
    #[serde(default = "default_output_dir")]
    pub output_dir: String,
    #[serde(default = "default_vtk_basename")]
    pub vtk_basename: String,
}

fn default_output_dir() -> String {
    "output".into()
}

fn default_vtk_basename() -> String {
    "windlab".into()
}

#[derive(Debug, Deserialize, Clone)]
pub struct GeometryConfig {
    /// Path to the STL mesh file.
    pub stl_path: String,
    #[serde(default = "default_padding")]
    pub padding: f64,
}

fn default_padding() -> f64 {
    0.2
}

impl SimConfig {
    pub fn load(path: impl AsRef<Path>) -> Result<Self, String> {
        let s = std::fs::read_to_string(path.as_ref()).map_err(|e| e.to_string())?;
        toml::from_str(&s).map_err(|e: toml::de::Error| e.to_string())
    }

    pub fn build_domain(&self) -> SoaDomain {
        let g = &self.grid;
        SoaDomain::new(g.nx, g.ny, g.nz, g.dx, g.dt)
    }

    pub fn lbm_params(&self) -> LbmParams {
        LbmParams {
            tau: self.physics.tau,
            body_force: self.physics.resolved_body_force(),
        }
    }

    pub fn apply_channel_walls_z(&self, types: &mut [NodeType]) {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let nz = self.grid.nz;
        for z in [0usize, nz - 1] {
            for y in 0..ny {
                for x in 0..nx {
                    let id = crate::lattice::index(nx, ny, x, y, z);
                    types[id] = NodeType::Solid;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn physics(tau: f64, body_force: [f64; 3]) -> PhysicsConfig {
        PhysicsConfig { tau, body_force }
    }

    #[test]
    fn body_force_array_when_no_magnitude() {
        let c = physics(0.9, [1.0, 2.0, 3.0]);
        assert_eq!(c.resolved_body_force(), [1.0, 2.0, 3.0]);
    }
}
