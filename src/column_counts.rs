use ndarray::{Array2, s};
use rayon::prelude::*;

/// Number of adjacent columns counted together, so each sequence is read in contiguous blocks.
pub const COL_BLOCK: usize = 64;

/// Number of times each byte value occurs in a column.
pub type ResidueCounts = [u32; 256];

/// Count the residues in each column of the alignment and summarise each column with `summarise`.
///
/// The alignment is stored row-major, so walking a single column touches one cache line per
/// sequence. Instead, blocks of `COL_BLOCK` adjacent columns are counted together by reading each
/// sequence's slice of the block, and blocks are processed in parallel. Summaries are returned in
/// column order.
pub fn map_column_counts<T, F>(alignment: &Array2<u8>, summarise: F) -> Vec<T>
where
    T: Send,
    F: Fn(&ResidueCounts) -> T + Sync,
{
    let ncols = alignment.ncols();
    let block_starts: Vec<usize> = (0..ncols).step_by(COL_BLOCK).collect();

    block_starts
        .par_iter()
        .map(|&start| {
            let end = (start + COL_BLOCK).min(ncols);
            let mut counts = vec![[0u32; 256]; end - start];
            for row in alignment.rows() {
                for (col_counts, &c) in counts.iter_mut().zip(row.slice(s![start..end])) {
                    col_counts[c as usize] += 1;
                }
            }
            counts.iter().map(&summarise).collect::<Vec<T>>()
        })
        .collect::<Vec<_>>()
        .into_iter()
        .flatten()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_map_column_counts() {
        let ncols = 2 * COL_BLOCK + 5;
        let alignment = Array2::from_shape_fn((3, ncols), |(r, c)| if r == 0 && c % 2 == 0 { b'-' } else { b'L' });
        let summaries = map_column_counts(&alignment, |counts| (counts[b'-' as usize], counts[b'L' as usize]));

        assert_eq!(summaries.len(), ncols);
        for (c, &s) in summaries.iter().enumerate() {
            let expected_gaps = if c % 2 == 0 { 1 } else { 0 };
            assert_eq!(s, (expected_gaps, 3 - expected_gaps));
        }
    }
}
