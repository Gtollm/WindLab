//! VTK ImageData (`.vti`) export via the `vtkio` crate

use std::path::Path;

use crate::geometry::coords::nm1;
use crate::geometry::stl::Bounds;
use crate::grid::cell::NodeType;
use crate::grid::SoaDomain;
use crate::lattice::CS2;

use vtkio::model::*;

pub fn write_vti_velocity(
    path: impl AsRef<Path>,
    domain: &SoaDomain,
    world_frame: Option<&Bounds>,
) -> std::io::Result<()> {
    let nx = domain.nx;
    let ny = domain.ny;
    let nz = domain.nz;
    let n = domain.ncells();

    let mut solid = vec![0.0_f32; n];
    let mut fluid = vec![0.0_f32; n];
    let mut rho = vec![0.0_f32; n];
    let mut pressure = vec![0.0_f32; n];
    let mut vel = Vec::<f32>::with_capacity(n * 3);

    let cs2 = CS2 as f32;
    let rho0 = 1.0_f32;

    for j in 0..n {
        if matches!(domain.node_type[j], NodeType::Solid) {
            solid[j] = 1.0;
            fluid[j] = 0.0;
            pressure[j] = 0.0;
            vel.push(0.0);
            vel.push(0.0);
            vel.push(0.0);
        } else {
            solid[j] = 0.0;
            fluid[j] = 1.0;
            let r = domain.rho[j] as f32;
            rho[j] = r;
            pressure[j] = cs2 * (r - rho0);
            vel.push(domain.ux[j] as f32);
            vel.push(domain.uy[j] as f32);
            vel.push(domain.uz[j] as f32);
        }
    }

    let (origin, spacing) = match world_frame {
        Some(b) => {
            let e = b.max - b.min;
            (
                [b.min.x as f32, b.min.y as f32, b.min.z as f32],
                [
                    (e.x / nm1(nx)) as f32,
                    (e.y / nm1(ny)) as f32,
                    (e.z / nm1(nz)) as f32,
                ],
            )
        }
        None => (
            [0.0; 3],
            [domain.dx as f32, domain.dx as f32, domain.dx as f32],
        ),
    };

    let extent = Extent::Ranges([
        0..=(nx as i32 - 1),
        0..=(ny as i32 - 1),
        0..=(nz as i32 - 1),
    ]);

    let point_data = Attributes {
        point: vec![
            Attribute::scalars("solid", 1).with_data(IOBuffer::new(solid)),
            Attribute::scalars("fluid", 1).with_data(IOBuffer::new(fluid)),
            Attribute::scalars("rho", 1).with_data(IOBuffer::new(rho)),
            Attribute::scalars("pressure", 1).with_data(IOBuffer::new(pressure)),
            Attribute::vectors("velocity").with_data(IOBuffer::new(vel)),
        ],
        cell: vec![],
    };

    let vtk = Vtk {
        version: Version::new((1, 0)),
        byte_order: ByteOrder::LittleEndian,
        title: String::from("WindLab LBM Simulation"),
        file_path: None,
        data: DataSet::ImageData {
            extent: extent.clone(),
            origin,
            spacing,
            meta: None,
            pieces: vec![Piece::Inline(Box::new(ImageDataPiece {
                extent,
                data: point_data,
            }))],
        },
    };

    vtk.export(path.as_ref())
        .map_err(|e| std::io::Error::other(e.to_string()))
}
