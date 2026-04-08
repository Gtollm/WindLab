//! VTK PolyData (`.vtp`) export of STL triangle surfaces

use std::path::Path;

use crate::geometry::stl::Tri;

use vtkio::model::*;

pub fn write_vtp_stl_surface(path: impl AsRef<Path>, tris: &[Tri]) -> std::io::Result<()> {
    let npolys = tris.len();
    let npts = npolys * 3;

    let mut points = Vec::with_capacity(npts * 3);
    for t in tris {
        for p in [t.a, t.b, t.c] {
            points.push(p.x as f32);
            points.push(p.y as f32);
            points.push(p.z as f32);
        }
    }

    let mut connectivity = Vec::with_capacity(npolys * 3);
    let mut offsets = Vec::with_capacity(npolys);
    for i in 0..npolys {
        let base = (i * 3) as u64;
        connectivity.push(base);
        connectivity.push(base + 1);
        connectivity.push(base + 2);
        offsets.push((i as u64 + 1) * 3);
    }

    let vtk = Vtk {
        version: Version::new((1, 0)),
        byte_order: ByteOrder::LittleEndian,
        title: String::from("WindLab STL Surface"),
        file_path: None,
        data: DataSet::inline(PolyDataPiece {
            points: IOBuffer::new(points),
            verts: None,
            lines: None,
            polys: Some(VertexNumbers::XML {
                connectivity,
                offsets,
            }),
            strips: None,
            data: Attributes::new(),
        }),
    };

    vtk.export(path.as_ref())
        .map_err(|e| std::io::Error::other(e.to_string()))
}
