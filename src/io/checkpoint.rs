//! Binary checkpoint save/restore for the SoA domain

use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

use crate::grid::cell::NodeType;
use crate::grid::SoaDomain;
use crate::lattice::Q;

const MAGIC: u32 = 0x574C_424D;

pub fn save_checkpoint(path: impl AsRef<Path>, domain: &SoaDomain) -> std::io::Result<()> {
    let mut f = File::create(path.as_ref())?;
    let n = domain.ncells();

    f.write_all(&MAGIC.to_le_bytes())?;
    f.write_all(&(domain.nx as u64).to_le_bytes())?;
    f.write_all(&(domain.ny as u64).to_le_bytes())?;
    f.write_all(&(domain.nz as u64).to_le_bytes())?;
    f.write_all(&domain.dx.to_le_bytes())?;
    f.write_all(&domain.dt.to_le_bytes())?;

    for j in 0..n {
        let tag: u8 = match domain.node_type[j] {
            NodeType::Fluid => 0,
            NodeType::Solid => 1,
            NodeType::Inlet(_) => 2,
            NodeType::Outlet(_) => 3,
        };
        f.write_all(&[tag])?;
    }

    for i in 0..Q {
        for j in 0..n {
            f.write_all(&domain.f[i][j].to_le_bytes())?;
        }
    }

    for j in 0..n {
        f.write_all(&domain.rho[j].to_le_bytes())?;
        f.write_all(&domain.ux[j].to_le_bytes())?;
        f.write_all(&domain.uy[j].to_le_bytes())?;
        f.write_all(&domain.uz[j].to_le_bytes())?;
    }

    Ok(())
}

pub fn load_checkpoint(path: impl AsRef<Path>, domain: &mut SoaDomain) -> std::io::Result<()> {
    let mut f = File::open(path.as_ref())?;
    let n = domain.ncells();

    let mut buf4 = [0u8; 4];
    f.read_exact(&mut buf4)?;
    if u32::from_le_bytes(buf4) != MAGIC {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "bad checkpoint magic",
        ));
    }

    let mut buf8 = [0u8; 8];
    f.read_exact(&mut buf8)?;
    let nx = u64::from_le_bytes(buf8) as usize;
    f.read_exact(&mut buf8)?;
    let ny = u64::from_le_bytes(buf8) as usize;
    f.read_exact(&mut buf8)?;
    let nz = u64::from_le_bytes(buf8) as usize;

    if nx != domain.nx || ny != domain.ny || nz != domain.nz {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "checkpoint grid mismatch: expected {} {} {}, found {} {} {}",
                domain.nx, domain.ny, domain.nz, nx, ny, nz
            ),
        ));
    }

    f.read_exact(&mut buf8)?;
    f.read_exact(&mut buf8)?;

    let mut tag = [0u8; 1];
    for j in 0..n {
        f.read_exact(&mut tag)?;
        domain.node_type[j] = match tag[0] {
            0 => NodeType::Fluid,
            1 => NodeType::Solid,
            _ => NodeType::Fluid,
        };
    }

    for i in 0..Q {
        for j in 0..n {
            f.read_exact(&mut buf8)?;
            domain.f[i][j] = f64::from_le_bytes(buf8);
        }
    }

    for j in 0..n {
        f.read_exact(&mut buf8)?;
        domain.rho[j] = f64::from_le_bytes(buf8);
        f.read_exact(&mut buf8)?;
        domain.ux[j] = f64::from_le_bytes(buf8);
        f.read_exact(&mut buf8)?;
        domain.uy[j] = f64::from_le_bytes(buf8);
        f.read_exact(&mut buf8)?;
        domain.uz[j] = f64::from_le_bytes(buf8);
    }

    for i in 0..Q {
        domain.f_tmp[i].copy_from_slice(&domain.f[i]);
    }

    Ok(())
}
