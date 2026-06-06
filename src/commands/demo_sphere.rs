use wind_lab::boundary::zou_he::{tag_x0_inlet, tag_xmax_outlet};
use wind_lab::core::solver::{step_soa, LbmParams};
use wind_lab::grid::NodeType;
use wind_lab::grid::SoaDomain;
use wind_lab::lattice::index;
use wind_lab::visualization::rerun_viz;

use super::utils::progress_bar;

pub fn run_demo_sphere(
    re: f64,
    diameter: usize,
    steps: usize,
    tau: f64,
    no_progress: bool,
    use_rerun: bool,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let nu = (tau - 0.5) / 3.0;
    let u_inlet = re * nu / diameter as f64;

    if u_inlet > 0.1 {
        eprintln!("Warning: u_inlet={u_inlet:.4} > 0.1 (Ma limit). Increase diameter or lower Re.");
    }

    let d = diameter;
    let nx = 6 * d;
    let ny = 8 * d;
    let nz = 8 * d;

    let mut domain = SoaDomain::new(nx, ny, nz, 1.0, 1.0);
    let mut types = vec![NodeType::Fluid; domain.ncells()];

    let cx = (1.5 * d as f64) as i32;
    let cy = (ny / 2) as i32;
    let cz = (nz / 2) as i32;
    let r2 = (d as f64 / 2.0).powi(2);
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let dx = x as i32 - cx;
                let dy = y as i32 - cy;
                let dz = z as i32 - cz;
                if (dx * dx + dy * dy + dz * dz) as f64 <= r2 {
                    types[index(nx, ny, x, y, z)] = NodeType::Solid;
                }
            }
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

    let rec = if use_rerun {
        let r = rerun_viz::spawn_viewer("windlab-sphere")?;
        rerun_viz::log_geometry(&r, &domain)?;
        Some(r)
    } else {
        None
    };

    let every = (steps / 40).max(1);
    let pb = progress_bar(steps, no_progress);
    for step in 0..steps {
        step_soa(&mut domain, &params);
        if let Some(p) = &pb {
            p.inc(1);
        }
        if let Some(r) = &rec {
            if (step + 1) % every == 0 {
                rerun_viz::log_drag(r, &domain, step + 1)?;
            }
        }
    }
    if let Some(p) = &pb {
        p.finish_with_message("sphere done");
    }

    let [fx, _, _] = wind_lab::physics::drag_force(&domain);
    let a_ref = wind_lab::physics::projected_area_yz(&domain).max(1) as f64;
    let cd_sim = wind_lab::physics::drag_coefficient(fx, 1.0, u_inlet, a_ref);
    let cd_ref = cd_empirical(re);

    let stokes_fx = 3.0 * std::f64::consts::PI * nu * d as f64 * u_inlet;
    let front_id = domain.idx((cx - 1) as usize, cy as usize, cz as usize);
    let back_id = domain.idx((cx + 1) as usize, cy as usize, cz as usize);
    let rho_front = domain.rho[front_id];
    let rho_back = domain.rho[back_id];

    println!("\n=== Sphere validation (Re={re:.1}) ===");
    println!("  tau        = {tau:.3}");
    println!("  u_inlet    = {u_inlet:.5}  (lattice units)");
    println!("  diameter   = {d} cells");
    println!("  grid       = {nx}x{ny}x{nz}");
    println!(
        "  Fx (LBM)   = {fx:.6}  Fx (Stokes) = {stokes_fx:.6}  ratio = {:.3}",
        fx / stokes_fx
    );
    println!(
        "  rho_front  = {rho_front:.6}  rho_back = {rho_back:.6}  d_rho = {:.6}",
        rho_front - rho_back
    );
    println!(
        "  a_ref      = {a_ref:.1}  (expected ~ {:.1})",
        std::f64::consts::PI * (d as f64 / 2.0).powi(2)
    );
    println!("  Cd (LBM)   = {cd_sim:.4}");
    println!("  Cd (ref)   = {cd_ref:.4}  [Oseen/empirical]");
    println!(
        "  error      = {:.1}%",
        (cd_sim - cd_ref).abs() / cd_ref * 100.0
    );

    Ok(())
}

fn cd_empirical(re: f64) -> f64 {
    if re < 1e-6 {
        return f64::INFINITY;
    }
    if re <= 1.0 {
        24.0 / re * (1.0 + 3.0 / 16.0 * re)
    } else if re <= 800.0 {
        24.0 / re * (1.0 + 0.15 * re.powf(0.687))
    } else {
        0.44
    }
}
