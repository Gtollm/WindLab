use wind_lab::core::solver::{step_soa, LbmParams};
use wind_lab::grid::cell::NodeType;
use wind_lab::grid::SoaDomain;
use wind_lab::io::write_vti_velocity;
use wind_lab::lattice::index;
use wind_lab::visualization::rerun_viz;

use super::utils::{parse_z_indices, progress_bar};

pub fn run_demo(
    n: usize,
    steps: usize,
    no_progress: bool,
    use_rerun: bool,
    viz_every: Option<usize>,
    slice_z_spec: Option<String>,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let mut domain = SoaDomain::new(n, n.max(4), n, 1.0, 1.0);
    let mut types = vec![NodeType::Fluid; domain.ncells()];
    let (nx, ny, nz) = (domain.nx, domain.ny, domain.nz);

    for z in [0usize, nz - 1] {
        for y in 0..ny {
            for x in 0..nx {
                types[index(nx, ny, x, y, z)] = NodeType::Solid;
            }
        }
    }

    domain.copy_node_types_from(&types);
    domain.init_uniform(1.0, 0.0, 0.0, 0.0);

    let params = LbmParams {
        tau: 0.9,
        body_force: [1e-5, 0.0, 0.0],
        boundary: Default::default(),
    };

    let rec = if use_rerun {
        let r = rerun_viz::spawn_viewer("windlab-demo")?;
        rerun_viz::log_geometry(&r, &domain)?;
        Some(r)
    } else {
        None
    };

    let z_indices = parse_z_indices(slice_z_spec.as_deref(), nz)?;
    let every = viz_every.unwrap_or((steps / 20).max(1));
    let pb = progress_bar(steps, no_progress);
    for step in 0..steps {
        step_soa(&mut domain, &params);
        if let Some(p) = &pb { p.inc(1); }
        if let Some(r) = &rec {
            if (step + 1) % every == 0 {
                rerun_viz::log_velocity_points(r, &domain, step + 1)?;
                rerun_viz::log_velocity_slice(r, &domain, step + 1, &z_indices)?;
                rerun_viz::log_drag(r, &domain, step + 1)?;
            }
        }
    }
    if let Some(p) = &pb { p.finish_with_message("demo finished"); }

    std::fs::create_dir_all("output")?;
    write_vti_velocity("output/demo_channel.vti", &domain, None)?;
    tracing::info!("Wrote output/demo_channel.vti ({steps} steps).");
    Ok(())
}
