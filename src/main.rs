use std::path::PathBuf;

use clap::{Parser, Subcommand};
use tracing::Level;
use tracing_subscriber::FmtSubscriber;

mod commands;

use commands::{
    demo_channel::run_demo, demo_sphere::run_demo_sphere, eval_stl::run_eval_stl,
    run::run_simulation,
};

fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let subscriber = FmtSubscriber::builder()
        .with_max_level(Level::INFO)
        .finish();
    let _ = tracing::subscriber::set_global_default(subscriber);

    Cli::parse().command.run()
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
        #[arg(long)]
        rerun: bool,
        #[arg(long)]
        slice_z: Option<String>,
    },
    DemoChannel {
        #[arg(long, default_value_t = 32)]
        n: usize,
        #[arg(long, default_value_t = 2000)]
        steps: usize,
        #[arg(long)]
        no_progress: bool,
        #[arg(long)]
        rerun: bool,
        #[arg(long)]
        viz_every: Option<usize>,
        #[arg(long)]
        slice_z: Option<String>,
    },
    EvalStl {
        #[arg(long)]
        stl: PathBuf,
        #[arg(long, default_value_t = 50.0)]
        re: f64,
        #[arg(long, default_value_t = 0.6)]
        tau: f64,
        #[arg(long, default_value_t = 30)]
        cpd: usize,
        #[arg(long, default_value_t = 8000)]
        steps: usize,
        #[arg(long)]
        no_progress: bool,
        #[arg(long)]
        quiet: bool,
    },
    DemoSphere {
        #[arg(long, default_value_t = 100.0)]
        re: f64,
        #[arg(long, default_value_t = 20)]
        diameter: usize,
        #[arg(long, default_value_t = 3000)]
        steps: usize,
        #[arg(long, default_value_t = 0.6)]
        tau: f64,
        #[arg(long)]
        no_progress: bool,
        #[arg(long)]
        rerun: bool,
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
                slice_z,
            } => run_simulation(&config, async_io, no_progress, rerun, slice_z),
            Commands::DemoChannel {
                n,
                steps,
                no_progress,
                rerun,
                viz_every,
                slice_z,
            } => run_demo(n, steps, no_progress, rerun, viz_every, slice_z),
            Commands::EvalStl {
                stl,
                re,
                tau,
                cpd,
                steps,
                no_progress,
                quiet,
            } => run_eval_stl(stl, re, tau, cpd, steps, no_progress, quiet),
            Commands::DemoSphere {
                re,
                diameter,
                steps,
                tau,
                no_progress,
                rerun,
            } => run_demo_sphere(re, diameter, steps, tau, no_progress, rerun),
        }
    }
}
