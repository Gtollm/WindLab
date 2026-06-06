use std::path::PathBuf;

use indicatif::{ProgressBar, ProgressStyle};

pub fn progress_bar(total: usize, disabled: bool) -> Option<ProgressBar> {
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

pub fn vtk_path(output_dir: &str, basename: &str, step: usize) -> PathBuf {
    PathBuf::from(output_dir).join(format!("{basename}_{step:06}.vti"))
}

pub fn parse_z_indices(
    spec: Option<&str>,
    nz: usize,
) -> Result<Vec<usize>, Box<dyn std::error::Error + Send + Sync>> {
    let Some(s) = spec else {
        return Ok(vec![nz / 2]);
    };

    let mut set = std::collections::BTreeSet::new();
    for token in s.split(',') {
        let token = token.trim();
        if token.is_empty() {
            continue;
        }
        if let Some((a, b)) = token.split_once('-') {
            let lo: usize = a
                .trim()
                .parse()
                .map_err(|_| format!("invalid index \"{a}\""))?;
            let hi: usize = b
                .trim()
                .parse()
                .map_err(|_| format!("invalid index \"{b}\""))?;
            if lo > hi {
                return Err(format!("range \"{token}\": start > end").into());
            }
            for z in lo..=hi {
                set.insert(z.min(nz - 1));
            }
        } else {
            let z: usize = token
                .parse()
                .map_err(|_| format!("invalid index \"{token}\""))?;
            set.insert(z.min(nz - 1));
        }
    }

    if set.is_empty() {
        return Ok(vec![nz / 2]);
    }
    Ok(set.into_iter().collect())
}
