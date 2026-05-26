use std::path::PathBuf;

use wind_lab::boundary::zou_he::{tag_x0_inlet, tag_xmax_outlet};
use wind_lab::core::solver::{step_soa, LbmParams};
use wind_lab::geometry::stl::Bounds;
use wind_lab::geometry::{load_stl_triangles, voxelize_triangles};
use wind_lab::grid::cell::NodeType;
use wind_lab::grid::SoaDomain;

use super::utils::progress_bar;

pub fn run_eval_stl(
    stl_path: PathBuf,
    re: f64,
    tau: f64,
    cpd: usize,
    steps: usize,
    no_progress: bool,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let (tris, stl_bounds) =
        load_stl_triangles(&stl_path).map_err(|e| format!("STL load: {e}"))?;

    let extent = stl_bounds.extent();
    let d_phys = extent.y.max(extent.z);
    if d_phys < 1e-12 {
        return Err("STL has zero cross-stream extent".into());
    }
    let dx = d_phys / cpd as f64;

    let front_pad = 1.5 * d_phys;
    let back_pad  = 3.0 * d_phys;
    let side_pad  = 1.0 * d_phys;

    let cy = (stl_bounds.min.y + stl_bounds.max.y) * 0.5;
    let cz = (stl_bounds.min.z + stl_bounds.max.z) * 0.5;

    let dom_min = nalgebra::Vector3::new(
        stl_bounds.min.x - front_pad,
        cy - (extent.y * 0.5 + side_pad),
        cz - (extent.z * 0.5 + side_pad),
    );
    let dom_max = nalgebra::Vector3::new(
        stl_bounds.max.x + back_pad,
        cy + (extent.y * 0.5 + side_pad),
        cz + (extent.z * 0.5 + side_pad),
    );
    let dom_extent = dom_max - dom_min;

    let nx = (dom_extent.x / dx).round() as usize + 1;
    let ny = (dom_extent.y / dx).round() as usize + 1;
    let nz = (dom_extent.z / dx).round() as usize + 1;

    let nu = (tau - 0.5) / 3.0;
    let u_inlet = re * nu / cpd as f64;
    if u_inlet > 0.1 {
        eprintln!(
            "Warning: u_inlet={u_inlet:.4} > 0.1 (Ma limit). Lower re or raise cpd."
        );
    }

    let world_bounds = Bounds { min: dom_min, max: dom_max };
    let mut domain = SoaDomain::new(nx, ny, nz, 1.0, 1.0);
    let mut types = vec![NodeType::Fluid; domain.ncells()];

    let vox = voxelize_triangles(nx, ny, nz, &tris, &world_bounds);
    for i in 0..types.len() {
        if matches!(vox[i], NodeType::Solid) {
            types[i] = NodeType::Solid;
        }
    }

    let u_vec = nalgebra::Vector3::new(u_inlet, 0.0, 0.0);
    tag_x0_inlet(&mut types, nx, ny, nz, u_vec);
    tag_xmax_outlet(&mut types, nx, ny, nz, 1.0);

    domain.copy_node_types_from(&types);
    domain.init_uniform(1.0, u_inlet, 0.0, 0.0);

    let params = LbmParams {
        tau,
        body_force: [0.0; 3],
        boundary: Default::default(),
    };

    let solid_count = types.iter().filter(|t| matches!(t, NodeType::Solid)).count();
    eprintln!("grid={nx}x{ny}x{nz}  solid={solid_count}  u_in={u_inlet:.5}  tau={tau:.3}");

    let pb = progress_bar(steps, no_progress);
    for _step in 0..steps {
        step_soa(&mut domain, &params);
        if let Some(p) = &pb { p.inc(1); }
    }
    if let Some(p) = &pb { p.finish_with_message("eval done"); }

    let [fx, _, _] = wind_lab::physics::drag_force(&domain);
    let a_ref = wind_lab::physics::projected_area_yz(&domain).max(1) as f64;
    let cd = wind_lab::physics::drag_coefficient(fx, 1.0, u_inlet, a_ref);

    eprintln!("Fx={fx:.6}  a_ref={a_ref:.1}");
    println!("Cd={cd:.6}");

    Ok(())
}
