//! STL triangle mesh loading

use std::fs::File;
use std::path::Path;

use nalgebra::Vector3;
use stl_io::read_stl;

#[derive(Clone, Debug)]
pub struct Tri {
    pub a: Vector3<f64>,
    pub b: Vector3<f64>,
    pub c: Vector3<f64>,
}

#[derive(Clone, Copy, Debug)]
pub struct Bounds {
    pub min: Vector3<f64>,
    pub max: Vector3<f64>,
}

impl Bounds {
    pub fn center(&self) -> Vector3<f64> {
        (self.min + self.max) * 0.5
    }

    pub fn extent(&self) -> Vector3<f64> {
        self.max - self.min
    }
}

pub fn load_stl_triangles(path: impl AsRef<Path>) -> Result<(Vec<Tri>, Bounds), String> {
    let mut f = File::open(path.as_ref()).map_err(|e| e.to_string())?;
    let mesh = read_stl(&mut f).map_err(|e| e.to_string())?;

    let mut tris = Vec::with_capacity(mesh.faces.len());
    let mut min = Vector3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut max = Vector3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);

    for face in &mesh.faces {
        let a = vertex_to_vec(&mesh.vertices[face.vertices[0]]);
        let b = vertex_to_vec(&mesh.vertices[face.vertices[1]]);
        let c = vertex_to_vec(&mesh.vertices[face.vertices[2]]);

        for p in [a, b, c] {
            min = min.inf(&p);
            max = max.sup(&p);
        }
        tris.push(Tri { a, b, c });
    }

    if tris.is_empty() {
        return Err("STL file contains no triangles".into());
    }

    Ok((tris, Bounds { min, max }))
}

fn vertex_to_vec(v: &stl_io::Vector<f32>) -> Vector3<f64> {
    Vector3::new(v[0] as f64, v[1] as f64, v[2] as f64)
}

pub fn pad_bounds(b: &mut Bounds, eps: f64) {
    let e = b.extent();
    if e.x < eps {
        b.min.x -= eps;
        b.max.x += eps;
    }
    if e.y < eps {
        b.min.y -= eps;
        b.max.y += eps;
    }
    if e.z < eps {
        b.min.z -= eps;
        b.max.z += eps;
    }
}

pub fn expand_bounds_relative(b: &mut Bounds, fraction: f64) {
    if fraction <= 0.0 {
        return;
    }
    let delta = b.extent() * fraction;
    b.min -= delta;
    b.max += delta;
}
