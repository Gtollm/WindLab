//! WindLab CLI - run LBM simulations from TOML config with optional VTK export

use std::path::PathBuf;

use clap::{Parser, Subcommand};
use indicatif::{ProgressBar, ProgressStyle};
use tracing::{info, Level};
use tracing_subscriber::FmtSubscriber;

use wind_lab::config::SimConfig;
use wind_lab::core::solver::{step_soa, LbmParams};
use wind_lab::geometry::stl::{expand_bounds_relative, pad_bounds, Bounds};
use wind_lab::geometry::{load_stl_triangles, voxelize_triangles};
use wind_lab::grid::cell::NodeType;
use wind_lab::grid::SoaDomain;
use wind_lab::io::{write_vti_velocity, write_vtp_stl_surface};
use wind_lab::visualization::rerun_viz::{self, RecordingStream};

fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    init_logging();
    let cli = Cli::parse();
    cli.command.run()
}

fn init_logging() {
    let subscriber = FmtSubscriber::builder()
        .with_max_level(Level::INFO)
        .finish();
    let _ = tracing::subscriber::set_global_default(subscriber);
}

#[derive(Parser)]
#[command(name = "windlab", version, about = "WindLab D3Q19 LBM CFD solver")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Run {
        #[arg(short, long)]
        config: PathBuf,
        #[arg(long, default_value_t = true)]
        async_io: bool,
        #[arg(long)]
        no_progress: bool,
        /// Stream live visualization to Rerun viewer (requires `rerun` CLI in PATH)
        #[arg(long)]
        rerun: bool,
    },
    DemoChannel {
        #[arg(long, default_value_t = 32)]
        n: usize,
        #[arg(long, default_value_t = 2000)]
        steps: usize,
        #[arg(long)]
        no_progress: bool,
        /// Stream live visualization to Rerun viewer (requires `rerun` CLI in PATH)
        #[arg(long)]
        rerun: bool,
        /// Log a Rerun frame every N steps (default: steps/20)
        #[arg(long)]
        viz_every: Option<usize>,
    },
}

impl Commands {
    fn run(self) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        match self {
            Commands::Run {
                config,
                async_io,
                no_progress,
                rerun,
            } => run_simulation(&config, async_io, no_progress, rerun),
            Commands::DemoChannel {
                n,
                steps,
                no_progress,
                rerun,
                viz_every,
            } => run_demo(n, steps, no_progress, rerun, viz_every),
        }
    }
}

fn run_simulation(
    config_path: &PathBuf,
    async_io: bool,
    no_progress: bool,
    use_rerun: bool,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let cfg = SimConfig::load(config_path)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

    let mut domain = cfg.build_domain();
    let params = cfg.lbm_params();

    let mut types = domain.node_type.clone();
    cfg.apply_channel_walls_z(&mut types);

    let (world_frame, stl_tris) = load_and_voxelize(&cfg, &mut types)?;

    domain.copy_node_types_from(&types);
    domain.init_uniform(1.0, 0.0, 0.0, 0.0);

    std::fs::create_dir_all(&cfg.io.output_dir)?;

    let rec = if use_rerun {
        let r = rerun_viz::spawn_viewer("windlab")?;
        if !stl_tris.is_empty() {
            if let Some(bounds) = &world_frame {
                rerun_viz::log_stl_mesh(&r, &stl_tris, bounds, cfg.grid.nx, cfg.grid.ny, cfg.grid.nz)?;
            }
        } else {
            rerun_viz::log_geometry(&r, &domain)?;
        }
        Some(r)
    } else {
        None
    };

    if async_io {
        run_async(&cfg, &params, &domain, world_frame, no_progress, rec)?;
    } else {
        run_sync(&cfg, &params, &domain, world_frame, no_progress, rec)?;
    }

    info!("Run finished.");
    Ok(())
}

fn load_and_voxelize(
    cfg: &SimConfig,
    types: &mut [NodeType],
) -> Result<(Option<Bounds>, Vec<wind_lab::geometry::stl::Tri>), Box<dyn std::error::Error + Send + Sync>> {
    let Some(geo) = &cfg.geometry else {
        return Ok((None, vec![]));
    };

    let (tris, mut bounds) = load_stl_triangles(&geo.stl_path).map_err(std::io::Error::other)?;

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
                rerun_viz::log_velocity_slice(r, &domain, step + 1)?;
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
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    use std::sync::Arc;
    use tokio::sync::Mutex;

    let rt = tokio::runtime::Runtime::new()?;
    rt.block_on(async {
        let out_dir = cfg.io.output_dir.clone();
        let base = cfg.io.vtk_basename.clone();
        let every = cfg.run.vtk_every;
        let steps = cfg.run.steps;
        let pb = progress_bar(steps, no_progress);

        let domain = Arc::new(Mutex::new(domain.clone()));
        let worker = domain.clone();

        for step in 0..steps {
            {
                let mut d = worker.lock().await;
                step_soa(&mut d, params);
            }
            if let Some(p) = &pb {
                p.inc(1);
            }
            if every > 0 && (step + 1) % every == 0 {
                let path = vtk_path(&out_dir, &base, step + 1);
                let snap = {
                    let d = worker.lock().await;
                    d.clone()
                };
                let frame = world_frame;
                let rec_snap = rec.clone();
                tokio::task::spawn_blocking(move || {
                    let _ = write_vti_velocity(&path, &snap, frame.as_ref());
                    if let Some(r) = rec_snap {
                        let _ = rerun_viz::log_velocity_points(&r, &snap, step + 1);
                        let _ = rerun_viz::log_velocity_slice(&r, &snap, step + 1);
                    }
                })
                .await
                .ok();
            }
        }

        if let Some(p) = &pb {
            p.finish_with_message("run finished");
        }
    });

    Ok(())
}

fn run_demo(
    n: usize,
    steps: usize,
    no_progress: bool,
    use_rerun: bool,
    viz_every: Option<usize>,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    use wind_lab::lattice::index;

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
    };

    let rec = if use_rerun {
        let r = rerun_viz::spawn_viewer("windlab-demo")?;
        rerun_viz::log_geometry(&r, &domain)?;
        Some(r)
    } else {
        None
    };

    let every = viz_every.unwrap_or((steps / 20).max(1));
    let pb = progress_bar(steps, no_progress);
    for step in 0..steps {
        step_soa(&mut domain, &params);
        if let Some(p) = &pb {
            p.inc(1);
        }
        if let Some(r) = &rec {
            if (step + 1) % every == 0 {
                rerun_viz::log_velocity_points(r, &domain, step + 1)?;
                rerun_viz::log_velocity_slice(r, &domain, step + 1)?;
            }
        }
    }
    if let Some(p) = &pb {
        p.finish_with_message("demo finished");
    }

    std::fs::create_dir_all("output")?;
    write_vti_velocity("output/demo_channel.vti", &domain, None)?;
    info!("Wrote output/demo_channel.vti ({steps} steps).");
    Ok(())
}

fn progress_bar(total: usize, disabled: bool) -> Option<ProgressBar> {
    if disabled || total == 0 {
        return None;
    }
    let pb = ProgressBar::new(total as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] \
                 {human_pos}/{human_len} steps (ETA {eta}) {msg}",
            )
            .expect("progress template")
            .progress_chars("=>-"),
    );
    pb.set_message("LBM");
    Some(pb)
}

fn vtk_path(output_dir: &str, basename: &str, step: usize) -> PathBuf {
    PathBuf::from(output_dir).join(format!("{basename}_{step:06}.vti"))
}
