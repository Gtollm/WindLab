use std::path::PathBuf;

use tracing::info;

use wind_lab::boundary::zou_he::{tag_x0_inlet, tag_xmax_outlet};
use wind_lab::config::SimConfig;
use wind_lab::core::solver::{step_soa, LbmParams};
use wind_lab::geometry::stl::{expand_bounds_relative, pad_bounds, Bounds};
use wind_lab::geometry::{load_stl_triangles, rotate_tris, voxelize_triangles};
use wind_lab::grid::NodeType;
use wind_lab::grid::SoaDomain;
use wind_lab::io::{write_vti_velocity, write_vtp_stl_surface};
use wind_lab::visualization::rerun_viz::{self, RecordingStream};

use super::utils::{parse_z_indices, progress_bar, vtk_path};

pub fn run_simulation(
    config_path: &PathBuf,
    async_io: bool,
    no_progress: bool,
    use_rerun: bool,
    slice_z_spec: Option<String>,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let cfg = SimConfig::load(config_path)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

    let mut domain = cfg.build_domain();
    let params = cfg.lbm_params();

    let mut types = domain.node_type.clone();
    cfg.apply_channel_walls_z(&mut types);

    let (world_frame, stl_tris) = load_and_voxelize(&cfg, &mut types)?;

    let iu = cfg.physics.inlet_u;
    if iu[0].abs() + iu[1].abs() + iu[2].abs() > f64::EPSILON {
        let u_vec = nalgebra::Vector3::new(iu[0], iu[1], iu[2]);
        tag_x0_inlet(&mut types, cfg.grid.nx, cfg.grid.ny, cfg.grid.nz, u_vec);
        tag_xmax_outlet(&mut types, cfg.grid.nx, cfg.grid.ny, cfg.grid.nz, 1.0);
    }

    domain.copy_node_types_from(&types);
    domain.init_uniform(1.0, 0.0, 0.0, 0.0);

    std::fs::create_dir_all(&cfg.io.output_dir)?;

    let rec = if use_rerun {
        let r = rerun_viz::spawn_viewer("windlab")?;
        if !stl_tris.is_empty() {
            if let Some(bounds) = &world_frame {
                rerun_viz::log_stl_mesh(
                    &r,
                    &stl_tris,
                    bounds,
                    cfg.grid.nx,
                    cfg.grid.ny,
                    cfg.grid.nz,
                )?;
            }
        } else {
            rerun_viz::log_geometry(&r, &domain)?;
        }
        Some(r)
    } else {
        None
    };

    let z_indices = parse_z_indices(slice_z_spec.as_deref(), cfg.grid.nz)?;
    if async_io {
        run_async(
            &cfg,
            &params,
            &domain,
            world_frame,
            no_progress,
            rec,
            z_indices,
        )?;
    } else {
        run_sync(
            &cfg,
            &params,
            &domain,
            world_frame,
            no_progress,
            rec,
            z_indices,
        )?;
    }

    info!("Run finished.");
    Ok(())
}

fn load_and_voxelize(
    cfg: &SimConfig,
    types: &mut [NodeType],
) -> Result<
    (Option<Bounds>, Vec<wind_lab::geometry::stl::Tri>),
    Box<dyn std::error::Error + Send + Sync>,
> {
    let Some(geo) = &cfg.geometry else {
        return Ok((None, vec![]));
    };

    let (mut tris, mut bounds) =
        load_stl_triangles(&geo.stl_path).map_err(std::io::Error::other)?;

    let rot = geo.rotation_deg;
    if rot[0].abs() + rot[1].abs() + rot[2].abs() > f64::EPSILON {
        bounds = rotate_tris(&mut tris, rot);
    }

    pad_bounds(&mut bounds, 1e-6);
    expand_bounds_relative(&mut bounds, geo.padding);

    let vtp_path =
        PathBuf::from(&cfg.io.output_dir).join(format!("{}_geometry.vtp", cfg.io.vtk_basename));
    write_vtp_stl_surface(&vtp_path, &tris)?;
    info!("STL surface overlay: {}", vtp_path.display());

    let vox = voxelize_triangles(cfg.grid.nx, cfg.grid.ny, cfg.grid.nz, &tris, &bounds);
    for i in 0..types.len() {
        if matches!(vox[i], NodeType::Solid) {
            types[i] = NodeType::Solid;
        }
    }

    Ok((Some(bounds), tris))
}

fn run_sync(
    cfg: &SimConfig,
    params: &LbmParams,
    domain: &SoaDomain,
    world_frame: Option<Bounds>,
    no_progress: bool,
    rec: Option<RecordingStream>,
    z_indices: Vec<usize>,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let mut domain = domain.clone();
    let pb = progress_bar(cfg.run.steps, no_progress);

    for step in 0..cfg.run.steps {
        step_soa(&mut domain, params);
        if let Some(p) = &pb {
            p.inc(1);
        }
        if cfg.run.vtk_every > 0 && (step + 1) % cfg.run.vtk_every == 0 {
            let path = vtk_path(&cfg.io.output_dir, &cfg.io.vtk_basename, step + 1);
            write_vti_velocity(&path, &domain, world_frame.as_ref())?;
            if let Some(r) = &rec {
                rerun_viz::log_velocity_points(r, &domain, step + 1)?;
                rerun_viz::log_velocity_slice(r, &domain, step + 1, &z_indices)?;
                rerun_viz::log_drag(r, &domain, step + 1)?;
            }
        }
    }

    if let Some(p) = &pb {
        p.finish_with_message("run finished");
    }
    Ok(())
}

fn run_async(
    cfg: &SimConfig,
    params: &LbmParams,
    domain: &SoaDomain,
    world_frame: Option<Bounds>,
    no_progress: bool,
    rec: Option<RecordingStream>,
    z_indices: Vec<usize>,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let rt = tokio::runtime::Runtime::new()?;
    rt.block_on(async {
        let out_dir = cfg.io.output_dir.clone();
        let base = cfg.io.vtk_basename.clone();
        let every = cfg.run.vtk_every;
        let steps = cfg.run.steps;
        let pb = progress_bar(steps, no_progress);

        let mut domain = domain.clone();
        let mut write_handle: Option<tokio::task::JoinHandle<()>> = None;

        for step in 0..steps {
            step_soa(&mut domain, params);
            if let Some(p) = &pb {
                p.inc(1);
            }

            if every > 0 && (step + 1) % every == 0 {
                if let Some(h) = write_handle.take() {
                    h.await.ok();
                }

                let snap = domain.clone();
                let path = vtk_path(&out_dir, &base, step + 1);
                let frame = world_frame;
                let rec_snap = rec.clone();
                let z_snap = z_indices.clone();

                write_handle = Some(tokio::task::spawn_blocking(move || {
                    let _ = write_vti_velocity(&path, &snap, frame.as_ref());
                    if let Some(r) = rec_snap {
                        let _ = rerun_viz::log_velocity_points(&r, &snap, step + 1);
                        let _ = rerun_viz::log_velocity_slice(&r, &snap, step + 1, &z_snap);
                        let _ = rerun_viz::log_drag(&r, &snap, step + 1);
                    }
                }));
            }
        }

        if let Some(h) = write_handle {
            h.await.ok();
        }
        if let Some(p) = &pb {
            p.finish_with_message("run finished");
        }
    });

    Ok(())
}
