//! D3Q19 lattice model constants and utilities.

pub mod equilibrium;

pub const Q: usize = 19;

pub const CS2: f64 = 1.0 / 3.0;

pub const C: [[i32; 3]; Q] = [
    [0, 0, 0],
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
    [1, 1, 0],
    [-1, 1, 0],
    [1, -1, 0],
    [-1, -1, 0],
    [1, 0, 1],
    [-1, 0, 1],
    [1, 0, -1],
    [-1, 0, -1],
    [0, 1, 1],
    [0, -1, 1],
    [0, 1, -1],
    [0, -1, -1],
];

pub const W: [f64; Q] = [
    1.0 / 3.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

pub const OPPOSITE: [usize; Q] = [
    0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15,
];

#[inline]
pub fn index(nx: usize, ny: usize, x: usize, y: usize, z: usize) -> usize {
    x + y * nx + z * nx * ny
}

#[inline]
pub fn unravel_index(nx: usize, ny: usize, id: usize) -> (usize, usize, usize) {
    let z = id / (nx * ny);
    let rem = id % (nx * ny);
    let y = rem / nx;
    let x = rem % nx;
    (x, y, z)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn opposite_pairs_sum_to_zero() {
        for i in 0..Q {
            let j = OPPOSITE[i];
            assert_eq!(C[i][0] + C[j][0], 0);
            assert_eq!(C[i][1] + C[j][1], 0);
            assert_eq!(C[i][2] + C[j][2], 0);
            assert_eq!(OPPOSITE[j], i);
        }
    }

    #[test]
    fn index_roundtrip() {
        let (nx, ny, nz) = (10, 8, 6);
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let id = index(nx, ny, x, y, z);
                    assert_eq!(unravel_index(nx, ny, id), (x, y, z));
                }
            }
        }
    }
}
