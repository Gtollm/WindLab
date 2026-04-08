//! Periodic boundary helpers for lattice indices

#[inline]
pub fn wrap(i: i32, n: usize) -> usize {
    let n = n as i32;
    let mut r = i % n;
    if r < 0 {
        r += n;
    }
    r as usize
}

fn neighbor_axis(s: i32, n: usize, periodic: bool) -> Option<usize> {
    if periodic {
        Some(wrap(s, n))
    } else if s >= 0 && (s as usize) < n {
        Some(s as usize)
    } else {
        None
    }
}

pub struct NeighborQuery {
    pub sx: i32,
    pub sy: i32,
    pub sz: i32,
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub px: bool,
    pub py: bool,
    pub pz: bool,
}

#[inline]
pub fn periodic_neighbor(q: NeighborQuery) -> Option<(usize, usize, usize)> {
    let gx = neighbor_axis(q.sx, q.nx, q.px)?;
    let gy = neighbor_axis(q.sy, q.ny, q.py)?;
    let gz = neighbor_axis(q.sz, q.nz, q.pz)?;
    Some((gx, gy, gz))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wrap_positive() {
        assert_eq!(wrap(5, 3), 2);
    }

    #[test]
    fn wrap_negative() {
        assert_eq!(wrap(-1, 3), 2);
        assert_eq!(wrap(-4, 3), 2);
    }

    #[test]
    fn periodic_neighbor_all_wrapped() {
        assert_eq!(
            periodic_neighbor(NeighborQuery {
                sx: -1,
                sy: 5,
                sz: 3,
                nx: 4,
                ny: 4,
                nz: 4,
                px: true,
                py: true,
                pz: true,
            }),
            Some((3, 1, 3))
        );
    }

    #[test]
    fn periodic_neighbor_out_of_bounds() {
        assert!(
            periodic_neighbor(NeighborQuery {
                sx: -1,
                sy: 5,
                sz: 3,
                nx: 4,
                ny: 4,
                nz: 4,
                px: false,
                py: true,
                pz: true,
            })
            .is_none()
        );
    }
}
