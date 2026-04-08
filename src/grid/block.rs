//! Block decomposition for parallel iteration along the Z axis

pub fn z_slab_ranges(nz: usize, num_chunks: usize) -> Vec<(usize, usize)> {
    if num_chunks == 0 || nz == 0 {
        return vec![(0, nz)];
    }
    let mut out = Vec::with_capacity(num_chunks);
    let base = nz / num_chunks;
    let rem = nz % num_chunks;
    let mut z = 0;
    for i in 0..num_chunks {
        let h = base + if i < rem { 1 } else { 0 };
        let z1 = (z + h).min(nz);
        out.push((z, z1));
        z = z1;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slabs_cover_grid() {
        let r = z_slab_ranges(10, 3);
        assert_eq!(r.len(), 3);
        let mut prev = 0;
        for (a, b) in &r {
            assert_eq!(*a, prev);
            prev = *b;
        }
        assert_eq!(prev, 10);
    }

    #[test]
    fn single_chunk_returns_full_range() {
        assert_eq!(z_slab_ranges(20, 1), vec![(0, 20)]);
    }

    #[test]
    fn zero_chunks_returns_single_range() {
        assert_eq!(z_slab_ranges(10, 0), vec![(0, 10)]);
    }
}
